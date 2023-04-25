function [] = corsequence(H)

S = guihandles(H);
D = guidata(H);

% test if there are images
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
img = {D.files.image};
% test if there is a ROI
if ~isfield(D,'roi')
    if isfield(D,'outputfile')
        headlessstatus('ERROR: no ROI defined, check the input file');
        return
    end
    msgstr = {'This action requires a defined ROI,';'set the ROI in section 3'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
roi = D.roi;

% test if there is a basis
if ~isfield(D,'basis')
    if isfield(D,'outputfile')
        headlessstatus('ERROR: no basis defined, check the input file');
        return
    end
    msgstr = {'This action requires a defined basis,';'set the basis in section 5'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end
basis = D.basis;

% load configuration options
cgstr{1,1} = 'prepA';
cgstr{2,1} = 'prepB';
cgstr{3,1} = 'prepC';
cgstr{4,1} = 'final';
for k = 1:4
    if k <= 3
        CG(k).use = get(S.([cgstr{k} 'on']),'Value');
    else
        CG(k).use = 1;
    end
    CG(k).blur = str2double(get(S.([cgstr{k} 'blur']),'String'));
    CG(k).cglevel = str2double(get(S.([cgstr{k} 'level']),'String'));
    CG(k).improcess = get(S.([cgstr{k} 'improcess']),'Value');
    CG(k).basis = get(S.([cgstr{k} 'basis']),'Value');
    CG(k).convcrit = str2double(get(S.([cgstr{k} 'convcrit']),'String'));
    CG(k).maxit = str2double(get(S.([cgstr{k} 'maxit']),'String'));
    CG(k).gradient = get(S.([cgstr{k} 'gradient']),'Value');
    CG(k).tikhpar1 = str2double(get(S.([cgstr{k} 'tikhpar1']),'String'));
    CG(k).tikhpar2 = str2double(get(S.([cgstr{k} 'tikhpar2']),'String'));
    CG(k).tikhsteps = str2double(get(S.([cgstr{k} 'tikhsteps']),'String'));
end

% more options
dimensions = get(S.dicdimensions,'Value')+1; % i.e. 2 or 3 for 2D or 3D
relaxation = get(S.dicrelaxation,'Value')-1; % i.e. 0-4, 0 = none, 1 = constant, etc.
relbasis = get(S.dicrelbasis,'Value');
lagrange = str2double(get(S.diclagrange,'String'));
convparam = get(S.dicconvparam,'Value'); % 1 du, 2 b,3 dp,4 dr
bestit = get(S.dicbestit,'Value'); % 1 r, 2 du, 3 b,4 dp, 5 dr, 6 lastit
maxdiv =  str2double(get(S.dicmaxdiv,'String'));
inciguess = get(S.dicinciguess,'Value'); % 1=0,2=iguess,3=prev, 4=iguess+prev, 5=reuse

autosave = get(S.dicautosave,'Value')-1;
memsave  = get(S.dicmemsave,'Value')-1;
precision = get(S.dicprecision,'Value');
usegpu = get(S.dicusegpu,'Value')-1;

liveview = get(S.corliveview,'Value')-1;
basedir = D.gui.ini.basedir.value;



% headlessmode
if isfield(D,'outputfile')
    headlessmode = true;
else
    headlessmode = false;
end

% check if all cg bases are done
% ===================
done = 0;
for k = 1:4
    if CG(k).use && isfield(basis,'type') && ~isempty(basis(CG(k).basis).type)
        done = done + 1;
    elseif ~CG(k).use
        done = done + 1;
    end
end
if done ~= 4
    if isfield(D,'outputfile')
        headlessstatus('ERROR: not all basis sets defined, check the input file');
        return
    end
    msgstr = {'This action requires a defined basis for each cg step,';'set the basis in section 5'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end

% list of basis functions
basislst = {D.basis.name}';

% update status
% ===================
D.gui.stat = appendstatus(D.gui.stat,'rule');
D.gui.stat = appendstatus(D.gui.stat,'[7] Correlation sequence started');
D.gui.stat = appendstatus(D.gui.stat,'rule');

% update log file
if headlessmode
    headlessstatus('rule')
    headlessstatus('[7] Correlation sequence started');
    headlessstatus('rule')
end

% load previous correlation results
if isfield(D,'cor')
    cor = D.cor;
else
    cor.inc = 1;
end

% reset status pane
set(S.corstatus,'String',{''});

% ask for autosave name
% ===================
if autosave
    if ~isempty(D.gui.ini.autosavefile.value)
        % save file specified in defaults.ini
        autosavefile = D.gui.ini.autosavefile.value;
        % convert relative path to absolute
        if ~strcmp(autosavefile(1),'/') && ~strcmp(autosavefile(2),':')
            autosavefile = fullfile(basedir,autosavefile);
        end
        if headlessmode
            headlessstatus(sprintf('autosaving each increment to %s',autosavefile));
        end
    elseif headlessmode
        % headless mode, use autputfile as autosave file
        autosavefile = D.outputfile;
        headlessstatus(sprintf('autosaving each increment to %s',autosavefile));
    else
        % ask for autosavefile location
        [filename, pathname] = uiputfile(fullfile(basedir,'autosave.mat'), 'autosave file name');
        if isequal(filename,0) || isequal(pathname,0)
            return
        else
            autosavefile = fullfile(pathname,filename);
        end
        
        % store the basename
        D.gui.ini.basedir.value = pathname;
    end
    
    % clear the old file (store one field)
    % this keeps the file nice and slim, otherwise matlab will accumulate
    % unneeded data in the file.
    version = D.version;
    save(autosavefile,'-v7.3','version');

    % open the matfile for writing
    mat = matfile(autosavefile,'Writable',true);
    
    % save to file
    mat.savetype = 'autosave';
    fields = fieldnames(D);
    for k = 1:length(fields)
        bcwaitbar(H,k/length(fields),sprintf('preparing autosave file'));
        mat.(fields{k}) = D.(fields{k});
    end
    
    % close the matfile
    delete(mat)
end

% test if GPU is compatible
% ===================
if usegpu
    try
        D.gpu = gpuDevice;
        guidata(H,D);
        % clear gpu memory
        reset(D.gpu)
    catch
        usegpu = 0;
        set(S.dicusegpu,'Value',1);
        if isfield(D,'outputfile')
            headlessstatus('WARNING: no compatible GPU detected, switching back to the CPU');
            return
        end
        msgstr = {'No compatible GPU device is found';'reverting to CPU computations only'};
        msgdlgjn(msgstr,dlgposition(H));
    end
end

% test memory availability
% ===================
if usegpu
    free = D.gpu.FreeMemory;
else
    free = memfree;
end

% estimated size of L
minsps = 2^min([CG.cglevel]);
npx = (roi(2)-roi(1)) * (roi(4)-roi(3));
ndof = 0;
for k = 1:4
    if CG(k).use
        ndof = max([ndof basis(CG(k).basis).Nphi]);
    end
end
if precision == 1
    % single precision
    precisionstr = 'single';
    byte = 4;
elseif precision == 2
    % double precision
    precisionstr = 'double';
    byte = 8;
end    
if memsave == 0
    est = (1/minsps)*npx*ndof*3*byte;
elseif memsave == 1
    est = (1/minsps)*npx*ndof*byte;
elseif memsave == 2
    est = (1/minsps)*npx*byte;
end
% safety factor (for other matrices besides L)
est = est * 2;

% status update
stat = sprintf('Current free memory %g Mb',round(free*1e-6));
D.gui.stat = appendstatus(D.gui.stat,stat);
stat = sprintf('Estimated required memory %g Mb',round(est*1e-6));
D.gui.stat = appendstatus(D.gui.stat,stat);
if headlessmode
    headlessstatus(sprintf('Current free memory %g Mb',round(free*1e-6)));
    headlessstatus(sprintf('Estimated required memory %g Mb',round(est*1e-6)));
end

% update memory stats to info
info(1).str = 'mem free (Mb)';
info(1).val = round(free*1e-6);
info(2).str = 'mem est (Mb)';
info(2).val = round(est*1e-6);
D.gui.info = appendinfo(D.gui.info,info);

% warn if estimated memory is more than free
if est > free
    qstr{1,1} = 'possible memory overflow detected';
    qstr{3,1} = sprintf('free memory: %g Mb',round(free*1e-6));
    qstr{4,1} = sprintf('estimated required memory: %g Mb',round(est*1e-6));
    qstr{6,1} = 'Continue anyway?';
    if isfield(D,'outputfile')
        headlessstatus(['WARNING: ' qstr{1}]);
        headlessstatus(['WARNING: ' qstr{2}]);
        headlessstatus(['WARNING: ' qstr{3}]);
        return
    end
    choice = questdlgjn(qstr,'possible memory overflow','Yes','No','No',dlgposition(H));
    if ~strcmp(choice,'Yes')
        D.gui.stat = appendstatus(D.gui.stat,'Aborted');
        guidata(H,D)
        return
    end
end

% clear previous basis matrices
% ===================
phi12{1} = [];
phi12{2} = [];
phi12{3} = [];
phi12{4} = [];
phi34{1} = [];
phi34{2} = [];
phi34{3} = [];
phi34{4} = [];
Ndims = {};

% initiate a convergence check list
convergecheck = [];
dic.converged = 0;

% number of increments and coarse grain steps
Ninc = length(img)-1;
Ncg = 4;
if length(cor) < Ninc
    cor(Ninc).inc = Ninc;
elseif length(cor) > Ninc
    cor = cor(1:Ninc);
end

% if zero increments, stop
if Ninc == 0
    return
end

% loop over the increments
for inc = 1:Ninc
    % update the waitbar
    bcwaitbar(H,inc/(Ninc+1),sprintf('correlating increment (%d/%d)',inc,Ninc));
    
    % update the gui
    set(S.corslider,'Value',inc);
    set(S.corid,'String',num2str(inc));

    % initiate increment structure
    cor(inc).inc = inc;
    
    % test if increment is already done
    if ~isfield(cor,'done')
        cor(inc).done = 0;
    end
    if isempty(cor(inc).done)
        cor(inc).done = 0;
    end
    if cor(inc).done == 4;
        % the increment is done, do not compute it again
        stat = get(S.corstatus,'String');
        stat = [sprintf('increment (%d) skipped',inc) ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        % status update
        stat = sprintf('  Skipping increment %d (already done)',inc);
        D.gui.stat = appendstatus(D.gui.stat,stat);
        if headlessmode
            headlessstatus(stat);
        end
        % skip this increment
        continue
    else
        % start the increment fresh
        cor(inc).done = 0;
    end
    
    % load images
    f = img{1};
    g = img{inc+1};
    [n, m] = size(f);
    x = 1:m;
    y = 1:n;
    
    % Updated Lagrange
    % ===========
    
    if (lagrange > 0) && (inc > 1) && (cor(inc-1).done == Ncg);
        % previous increment has no data yet
        stat = get(S.corstatus,'String');
        stat = [sprintf('=====================') ; stat];
        stat = ['lagrange, computing g-tilde' ; stat];
        stat = [sprintf('=====================') ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        % status update
        stat = 'Lagrange, computing the back deformed image g-tilde';
        D.gui.stat = appendstatus(D.gui.stat,stat);
        if headlessmode
            headlessstatus(stat);
        end
        
        % get the displacement fields
        if isfield(cor(inc-1),'U1') && ~isempty(cor(inc-1).U1)
            U1 = cor(inc-1).U1;
            U2 = cor(inc-1).U2;
        else
            error('unexpected error, no data to create lagrange image g-tilde');
        end
        Usiz = size(U1);
        U3 = zeros(Usiz);
        U4 = zeros(Usiz);
        U5 = zeros(Usiz);
        U6 = zeros(Usiz);
        if isfield(cor(inc-1),'U3') && ~isempty(cor(inc-1).U3)
            U3 = cor(inc-1).U3;
        end
        if isfield(cor(inc-1),'U4') && ~isempty(cor(inc-1).U4)
            U4 = cor(inc-1).U4;
        end
        if isfield(cor(inc-1),'U5') && ~isempty(cor(inc-1).U5)
            U5 = cor(inc-1).U5;
        end
        if isfield(cor(inc-1),'U6') && ~isempty(cor(inc-1).U6)
            U6 = cor(inc-1).U6;
        end

        % get the new coordinates
        [X, Y] = meshgrid(x,y);
        Im = find(x >= roi(1) & x <= roi(2));
        In = find(y >= roi(3) & y <= roi(4));
        Iroi = (X >= roi(1) & X <= roi(2) & Y >= roi(3) & Y <= roi(4));
        xroi = x(Im);
        yroi = y(In);
        [Xroi, Yroi] = meshgrid(xroi,yroi);
        
        % test if the last step is coarse grained
        if CG(Ncg).cglevel > 0
            % get the old coordinates
            xroi1 = cor(inc-1).xroi;
            yroi1 = cor(inc-1).yroi;
            [Xroi1, Yroi1] = meshgrid(xroi1,yroi1);
            
            % interpolate the displacement fields
            U1 = interp2(Xroi1,Yroi1,U1,Xroi,Yroi,'linear');
            U2 = interp2(Xroi1,Yroi1,U2,Xroi,Yroi,'linear');
            U3 = interp2(Xroi1,Yroi1,U3,Xroi,Yroi,'linear');
            U4 = interp2(Xroi1,Yroi1,U4,Xroi,Yroi,'linear');
            U5 = interp2(Xroi1,Yroi1,U5,Xroi,Yroi,'linear');
            U6 = interp2(Xroi1,Yroi1,U6,Xroi,Yroi,'linear');
        end
        
        % compute g-tilde
        gt = interp2(X,Y,img{inc},X(In,Im)+U1,Y(In,Im)+U2,'spline',NaN);
        % correct the brightness
        gt = gt + U3 + U4.*gt + U5.*gt.^2 + U6.*gt.^3;
        
        % if any NaN's occured, replace those with pixels from f
        gt(isnan(gt)) = f(isnan(gt));
        
        % mixing total and updated images
        f(Iroi) = (1-lagrange)*f(Iroi) + lagrange*gt(:);
        f = reshape(f,n,m);
        
    elseif (lagrange > 0) && (inc > 1)
        % previous increment has no data yet
        stat = get(S.corstatus,'String');
        stat = ['Skipping inc, prev. inc. is not done (lagrange)' ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        % status update
        stat = 'Skipping increment, previous increment is not done (lagrange)';
        D.gui.stat = appendstatus(D.gui.stat,stat);
        if headlessmode
            headlessstatus(stat);
        end
        % skip this increment
        continue
        
    end
    
    % force many matrices to single
    if precision == 1
        f = single(f);
        g = single(g);
    end
    
    % scale the images (for better M balancing)
    if dimensions == 2 % 2D
        % center the grayvalues
        meanf = mean(f(:));
        f = f - meanf;
        g = g - meanf;

        % scale with the dynamic image range
        drange = max(f(:)) - min(f(:));
        f = f./drange;
        g = g./drange;
    else % 3D
        % no scaling nor centering
        drange = 1;
        meanf = 0;
    end
    
    % mask image
    maskimg = maskimage(D.mask,n,m);
    
    for icg = 1:Ncg
        
        % status update
        stat = sprintf('increment: %d/%d, coarse grain step: %d/%d',inc,Ninc,icg,Ncg);
        D.gui.stat = appendstatus(D.gui.stat,'smallrule');
        D.gui.stat = appendstatus(D.gui.stat,stat);
        D.gui.stat = appendstatus(D.gui.stat,'smallrule');
        if headlessmode
            headlessstatus('smallrule');
            headlessstatus(stat);
            headlessstatus('smallrule');
        end
        
        % load cg settings
        cg.icg = icg;
        cg.Ncg = Ncg;
        cg.inc = inc;
        cg.Ninc = Ninc;
        cg.f = f;
        cg.g = g;
        cg.x = x;
        cg.y = y;
        cg.maskimg = maskimg;
        
        cg.use = CG(icg).use;
        cg.blur = CG(icg).blur;
        cg.level = CG(icg).cglevel;
        cg.sps = 2^cg.level;
        cg.improcess = CG(icg).improcess;
        cg.basisname = basislst{CG(icg).basis};
        cg.basisval = CG(icg).basis;
        cg.convcrit = CG(icg).convcrit;
        cg.maxit = CG(icg).maxit;
        cg.gradient = CG(icg).gradient;
        cg.tikhsteps = CG(icg).tikhsteps;
        cg.tikhpar1 = CG(icg).tikhpar1;
        cg.tikhpar2 = CG(icg).tikhpar2;
        
        cg.convparam = convparam;
        cg.memsave = memsave;
        cg.precision = precision;
        cg.maxdiv = maxdiv;
        cg.bestit = bestit;
        cg.dimensions = dimensions;
        cg.roi = roi;
        cg.liveview = liveview;
        cg.headlessmode = headlessmode;
        if relbasis == 1
            cg.relbasisname = cg.basisname;
            cg.relbasisval = cg.basisval;
        else
            cg.relbasisname = basislst{relbasis-1};
            cg.relbasisval = relbasis-1;
        end
        
        if ~cg.use
            % status update
            stat = sprintf('Skipping coarse grain step %d (disabled in options)',icg);
            D.gui.stat = appendstatus(D.gui.stat,stat);
            if headlessmode
                headlessstatus(stat);
            end
            continue
        end
        
        % status update
        stat = get(S.corstatus,'String');
        stat = [sprintf('=====================') ; stat];
        stat = [sprintf('inc %d/%d, cg %d/%d',inc,Ninc,icg,Ncg) ; stat];
        stat = [sprintf('=====================') ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        % Image Processing
        % ===========

        % update status
        stat = get(S.corstatus,'String');
        stat = ['pre-correlation (image processing)' ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        % blur
        [cg.f, cg.g] = improcessblur(cg);
        
        % gradient transform
        [cg.f, cg.g] = improcessgradient(cg);
        
        % coarsegrain
        [cg.f, cg.x, cg.y] = coarsegrain(cg.f,cg.level);
        [cg.g]             = coarsegrain(cg.g,cg.level);
        [cg.maskimg]       = coarsegrain(cg.maskimg,cg.level,0.25);
        [cg.n, cg.m]       = size(cg.f);
        [X, Y]             = meshgrid(cg.x,cg.y);
        
        
        % basis
        % ===========

        % update status
        stat{1} = 'pre-correlation (computing basis)';
        set(S.corstatus,'String',stat);drawnow
        
        if isempty(phi12{icg})
            
            % create the displacement basis matrix
            % x=>1st direction
            % y=>2nd direction
            phi12{icg} = phibuild(cg.x,cg.y,roi,basis(cg.basisval),H,precisionstr);
            phi34{icg} = [];
            
            if dimensions == 3 % 'Quasi 3D'
                % z=>3rd direction
                Ndims{icg} = 3;
                phi34{icg} = phi12{icg};
            elseif dimensions == 2 %2D'
                % brightness=>3rd direction 
                % contrast=>4th direction
                Ndims{icg} = relaxation + 2;
                
                % create the relaxation basis matrix
                if Ndims{icg} >= 3
                    if cg.relbasisval == cg.basisval
                        phi34{icg} = phi12{icg};
                    else
                        phi34{icg} = phibuild(cg.x,cg.y,roi,basis(cg.relbasisval),H,precisionstr);
                    end
                end

            end
            
            % update info panel
            Iroi = find(X(:) >= roi(1) & X(:) <= roi(2) & Y(:) >= roi(3) & Y(:) <= roi(4));
            Npx = numel(Iroi);
            Ndof = size(phi12{icg},2);
            Ntot = sum(phi12{icg}(Iroi,:)~=0);
            info(1).str = sprintf('phi12 cg%d Npx',icg);
            info(1).val = Npx;
            info(2).str = sprintf('phi12 cg%d Ndof',icg);
            info(2).val = Ndof;
            info(3).str = sprintf('phi12 cg%d Npx/Ndof',icg);
            info(3).val = Npx/Ndof;
            info(4).str = sprintf('phi12 cg%d support',icg);
            info(4).val = mean(Ntot);
            if Ndims{icg} >= 3
                Ndof = size(phi34{icg},2);
                Ntot = sum(phi34{icg}(Iroi,:)~=0);
                info(5).str = sprintf('phi34 cg%d Npx',icg);
                info(5).val = Npx;
                info(6).str = sprintf('phi34 cg%d Ndof',icg);
                info(6).val = Ndof;
                info(7).str = sprintf('phi34 cg%d Npx/Ndof',icg);
                info(7).val = Npx/Ndof;
                info(8).str = sprintf('phi34 cg%d support',icg);
                info(8).val = mean(Ntot);
            end
            D.gui.info = appendinfo(D.gui.info,info);
        end
        
        % store the number of dimensions to use in correlate.m
        cg.Ndims = Ndims{icg};
        
        % number of dof (per direction)
        cg.Nphi([1,2]) = size(phi12{icg},2);
        if cg.Ndims == 3
            cg.Nphi(3) = size(phi34{icg},2);
        elseif cg.Ndims == 4
            cg.Nphi([3,4]) = size(phi34{icg},2);
        elseif cg.Ndims == 5
            cg.Nphi([3,4,5]) = size(phi34{icg},2);
        elseif cg.Ndims == 6
            cg.Nphi([3,4,5,6]) = size(phi34{icg},2);
        end

        % initial guess
        % ===========
        
        % update status
        stat{1} = 'pre-correlation (initial guess)';
        set(S.corstatus,'String',stat);drawnow
        
        % initial guess interpolation status string
        pstr{1} = '';
        pstr{2} = '';
        
        if inc == 1
            % first increment
            if cor(inc).done == 0
                % first cg step
                if inciguess == 1 % 'zero'
                    cg.p = zeros(sum(cg.Nphi),1);
                elseif inciguess == 2 % init. guess
                    [cg.p, pstr] = piguess(H,D,cg,phi12{icg},phi34{icg});
                elseif inciguess == 3 % prev. inc.
                    cg.p = zeros(sum(cg.Nphi),1);
                elseif inciguess == 4 % prev. inc. and init. guess
                    [cg.p, pstr] = piguess(H,D,cg,phi12{icg},phi34{icg});
                elseif inciguess == 5 % reuse
                    if isfield(cor,'U1') && ~isempty(cor(inc).U1)
                        cg.p = pprev(H,cor(inc),cg,phi12{icg},phi34{icg});
                    else
                        cg.p = zeros(sum(cg.Nphi),1);
                    end
                end
            else
                % consecutive cg steps
                cg.p = pprev(H,dic,cg,phi12{icg},phi34{icg});
            end
        elseif inc > 1
            % consecutive increments
            if cor(inc).done == 0
                % first cg step
                if inciguess == 1 % 'zero'
                    cg.p = zeros(sum(cg.Nphi),1);
                elseif inciguess == 2 % init. guess
                    [cg.p, pstr] = piguess(H,D,cg,phi12{icg},phi34{icg});
                elseif inciguess == 3 % prev. inc.
                    cg.p = pprev(H,cor(inc-1),cg,phi12{icg},phi34{icg});
                elseif inciguess == 4 % prev. inc. and init. guess
                    [pi, pstr] = piguess(H,D,cg,phi12{icg},phi34{icg});
                    pp = pprev(H,cor(inc-1),cg,phi12{icg},phi34{icg});
                    cg.p = 0.5*(pi + pp);
                elseif inciguess == 5 % reuse
                    if isfield(cor,'U1') && ~isempty(cor(inc).U1)
                        cg.p = pprev(H,cor(inc),cg,phi12{icg},phi34{icg});
                    else
                        cg.p = pprev(H,cor(inc-1),cg,phi12{icg},phi34{icg});
                    end
                end
            else
                % consecutive cg steps
                cg.p = pprev(H,dic,cg,phi12{icg},phi34{icg});
            end
        end
        
        % update status
        stat = stat(2:end);
        set(S.corstatus,'String',stat);drawnow
        
        % initial guess interpolation message
        if ~isempty(pstr{1})
            % status update
            stat = sprintf('iguess: too few initial guess points, interpolating (Ux,Uy) with a %s order polynomial',pstr{1});
            D.gui.stat = appendstatus(D.gui.stat,stat);
            if cg.headlessmode
                headlessstatus(stat);
            end
        end
        if ~isempty(pstr{2})
            % status update
            stat = sprintf('iguess: too few initial guess points, interpolating (Uz) with a %s order polynomial',pstr{2});
            D.gui.stat = appendstatus(D.gui.stat,stat);
            if cg.headlessmode
                headlessstatus(stat);
            end
        end

        if get(S.corstop,'Value')
            stat = get(S.corstatus,'String');
            stat = [sprintf('coarse grain step stopped') ; stat];
            set(S.corstatus,'String',stat);drawnow
            D.gui.stat = appendstatus(D.gui.stat,'! Coarse grain step stopped !');
            break
        end
        
        % correlate
        % ===========

        dic = correlate(H,S,cg,phi12{icg},phi34{icg});

%         assignin('base','cg',cg);
%         assignin('base','dic',dic);
%         assignin('base','cor',cor);
%         assignin('base','D',D);
%         assignin('base','phi12',phi12);
%         assignin('base','phi34',phi34);
        
        if get(S.corstop,'Value')
            stat = get(S.corstatus,'String');
            stat = [sprintf('coarse grain step stopped') ; stat];
            set(S.corstatus,'String',stat);drawnow
            D.gui.stat = appendstatus(D.gui.stat,'! Coarse grain step stopped !');
            break
        end
        
        % append status
        D.gui.stat = [dic.statstr ; D.gui.stat];

        % add one to the done counter (done == 4, means all cgs are done)
        cor(inc).done = icg;
        
        % store
        % ===========
        
        % step data
        cor(inc).cg(icg).p  = dic.p;
        cor(inc).cg(icg).ip = cg.p;
        cor(inc).cg(icg).itable = single(dic.itable);
        
        % Results, only stored to cor for the final step
        if icg == Ncg
            cor(inc).x = cg.x;
            cor(inc).y = cg.y;
            cor(inc).In = dic.In;
            cor(inc).Im = dic.Im;
            cor(inc).xroi = dic.xroi;
            cor(inc).yroi = dic.yroi;
            cor(inc).Imask = dic.Imask;
            cor(inc).Iunmask = dic.Iunmask;

            cor(inc).meanf = meanf;
            cor(inc).drange = drange;
            cor(inc).converged = dic.converged;
            
            cor(inc).froi = dic.froi;
            cor(inc).r = dic.r;
            
            cor(inc).U1 = dic.U1;
            cor(inc).U2 = dic.U2;
            if cg.Ndims >= 3
                cor(inc).U3 = dic.U3;
            end
            if cg.Ndims >= 4
                cor(inc).U4 = dic.U4;
            end
            if cg.Ndims >= 5
                cor(inc).U5 = dic.U5;
            end
            if cg.Ndims == 6
                cor(inc).U6 = dic.U6;
            end
        end

    end % end cg loop
    
    if get(S.corstop,'Value')
        stat = get(S.corstatus,'String');
        stat = [sprintf('increment stopped') ; stat];
        set(S.corstatus,'String',stat);drawnow
        D.gui.stat = appendstatus(D.gui.stat,'! Increment stopped !');
        break
    end
    
    % autosave the current increment
    if autosave
        
        % open the matfile for writing
        mat = matfile(autosavefile,'Writable',true);
        
        % save to file
        mat.cor(1,inc) = cor(inc);
        mat.gui = D.gui;
        
        % close the matfile
        delete(mat)

        stat = get(S.corstatus,'String');
        stat = [sprintf('increment (%d) saved',inc) ; stat];
        set(S.corstatus,'String',stat);drawnow
        if headlessmode
            headlessstatus(stat(1));
        end
        
        
    end
    
    % add to the convergence checklist
    
    if dic.converged ~= 1
        convergecheck = [convergecheck inc];
    end
    
end % end inc loop
bcwaitbar(H);

% store the correlation results to the D structure
D.cor = cor;

% assignin('base','cor',cor);
% assignin('base','D',D);
% assignin('base','phi12',phi12);
% assignin('base','phi34',phi34);

% create a convergecheck str
convchkstr = sprintf(',%d',convergecheck);
convchkstr = convchkstr(2:end);

% Show if any increment failed
msgstr{1,1} = 'Correlation done';
if ~isempty(convergecheck)
    if length(convergecheck) == 1
        msgstr{2,1} = sprintf('WARNING: increment %s did not converge properly and should be checked',convchkstr);
    else
        msgstr{2,1} = sprintf('WARNING: increments [%s] did not converge properly and should be checked',convchkstr);
    end
else
    msgstr{2,1} = 'all increments converged properly';
end
if ~isfield(D,'outputfile') && ~get(S.corstop,'Value')
    % only if not headless or if not stopped
    msgdlgjn(msgstr,dlgposition(H));
end


stat = get(S.corstatus,'String');
stat = [sprintf('\ncorrelation done\n') ; stat];
set(S.corstatus,'String',stat);drawnow

% status update
D.gui.stat = appendstatus(D.gui.stat,'[7] Correlation done');
if ~isempty(convergecheck)
    if length(convergecheck) == 1
        msgstr = sprintf('[7] WARNING: increment %s did not converge properly and should be checked',convchkstr);
    else
        msgstr = sprintf('[7] WARNING: increments [%s] did not converge properly and should be checked',convchkstr);
    end
    D.gui.stat = appendstatus(D.gui.stat,msgstr);
end
D.gui.stat = appendstatus(D.gui.stat,'rule');
D.gui.stat = appendstatus(D.gui.stat,'');
D.gui.stat = appendstatus(D.gui.stat,'');
if headlessmode
    headlessstatus('rule');
    headlessstatus('[7] Correlation done');
    headlessstatus('rule');
end


% update the application data
guidata(H,D);


