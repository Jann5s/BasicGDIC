function dic = correlate(H,S,cg,phi12,phi34)
D = guidata(H);

usegpu = get(S.dicusegpu,'Value')-1;

if isfield(D,'outputfile')
    headlessmode = true;
else
    headlessmode = false;
end

% Prepare DIC
% ========================
dic = [];

% image space
f = cg.f;
g = cg.g;
x = cg.x;
y = cg.y;
[X, Y] = meshgrid(x,y);

if usegpu
    f = gpuArray(f);
    g = gpuArray(g);
    phi12 = gpuArray(phi12);
    phi34 = gpuArray(phi34);
end

% region of interest
roi = cg.roi;
Im = find(x >= roi(1) & x <= roi(2));
In = find(y >= roi(3) & y <= roi(4));
Iroi = find(X(:) >= roi(1) & X(:) <= roi(2) & Y(:) >= roi(3) & Y(:) <= roi(4));
xroi = x(Im);
yroi = y(In);
mroi = length(xroi);
nroi = length(yroi);
froi = f(In,Im);
[Xroi,Yroi] = meshgrid(xroi,yroi);

% advanced options
cornquiver    = D.gui.ini.cornquiver.value;
cordofsupport = D.gui.ini.cordofsupport.value;
corsparseness = D.gui.ini.corsparseness.value;
corgradcrit   = D.gui.ini.corgradcrit.value;
cordivtol     = D.gui.ini.cordivtol.value;

% tikhonov parameters
cg.tikhpar = logspace(log10(cg.tikhpar1),log10(cg.tikhpar2),cg.tikhsteps);

% initial guess
p = cg.p;
Ndof = length(p);

% number of dimensions
Ndims = cg.Ndims;
Idims{1} = 1:cg.Nphi(1);
Idims{2} = Idims{1}(end) + (1:cg.Nphi(2));
if Ndims >= 3
    Idims{3} = Idims{2}(end) + (1:cg.Nphi(3));
end
if Ndims >= 4
    Idims{4} = Idims{3}(end) + (1:cg.Nphi(4));
end
if Ndims >= 5
    Idims{5} = Idims{4}(end) + (1:cg.Nphi(5));
end
if Ndims >= 6
    Idims{6} = Idims{5}(end) + (1:cg.Nphi(6));
end

Idofdim = []; % lists which dimension this dof is for
Idofdof = []; % lists which subdof in phi12 or phi23
for k = 1:length(cg.Nphi)
    Idofdim = [Idofdim k*ones(1,cg.Nphi(k))];
    Idofdof = [Idofdof 1:cg.Nphi(k)];
end

% image gradients
dx = mean(diff(xroi));
dy = mean(diff(yroi));
[dfdx, dfdy] = gradient(f(In,Im),dx,dy);
[dgdx, dgdy] = gradient(g,dx,dy);

% get total support
support = sum(abs(phi12(Iroi,:)),2);
if Ndims >= 3
    support = support + sum(abs(phi34(Iroi,:)),2);
end
support = reshape(support,nroi,mroi);

% get mask
Imask = find(cg.maskimg(In,Im) == 1 | support == 0);
Iunmask = find(cg.maskimg(In,Im) ~= 1 & support > 0);
Npx = length(Iunmask);

% Quiver
Nq = cornquiver;
Iqn = round(linspace(2,nroi-1,Nq));
Iqm = round(linspace(2,mroi-1,Nq));
[IQm, IQn] = meshgrid(Iqm,Iqn);
Iq = sub2ind([nroi mroi],IQn(:),IQm(:));
Iq = intersect(Iq,Iunmask);

% liveviewdata
lview.xlim = [x(1) x(end)];
lview.ylim = [y(1) y(end)];
lview.xroi = xroi;
lview.yroi = yroi;
lview.X = Xroi(Iq);
lview.Y = Yroi(Iq);
lview.r = froi;
lview.U1 = zeros(length(Iq),1);
lview.U2 = zeros(length(Iq),1);
lview.inc = cg.inc;
lview.Ninc = cg.Ninc;
lview.icg = cg.icg;
lview.it = 0;

% Check for ill-supported dof (as cons. of masking)
Ntot(Idims{1}) = sum(phi12(Iroi,:)~=0);
Nsup(Idims{1}) = sum(phi12(Iroi(Iunmask),:)~=0);
Ntot(Idims{2}) = Ntot(Idims{1});
Nsup(Idims{2}) = Nsup(Idims{1});
if Ndims >= 3
    Ntot(Idims{3}) = sum(phi34(Iroi,:)~=0);
    Nsup(Idims{3}) = sum(phi34(Iroi(Iunmask),:)~=0);
end
if Ndims >= 4
    Ntot(Idims{4}) = Ntot(Idims{3});
    Nsup(Idims{4}) = Nsup(Idims{3});
end
if Ndims >= 5
    Ntot(Idims{5}) = Ntot(Idims{3});
    Nsup(Idims{5}) = Nsup(Idims{3});
end
if Ndims == 6
    Ntot(Idims{6}) = Ntot(Idims{3});
    Nsup(Idims{6}) = Nsup(Idims{3});
end

% indices of the well supported shape functions
Isup = find(Nsup./Ntot > cordofsupport);
Iunsup = find(Nsup./Ntot <= cordofsupport);

% update status
stat = get(S.corstatus,'String');
stat = [sprintf('(%3d/%3d) free dof',length(Isup),sum(cg.Nphi)) ; stat];
set(S.corstatus,'String',stat);drawnow

if ~isempty(Iunsup)
    for k = 1:length(Iunsup);
        if k == 1
            unsupportedstr = num2str(Iunsup(k));
        else
            unsupportedstr = [unsupportedstr, ',' num2str(Iunsup(k))];
        end
        if mod(k,10) == 0
            unsupportedstr = [unsupportedstr ' '];
        end
    end
    str = sprintf('(%3d/%3d) free dof, unsupported: [%s]',length(Isup),sum(cg.Nphi),unsupportedstr);
    statstr = appendstatus({''},str);
    if cg.headlessmode
        headlessstatus(str);
    end
else
    str = sprintf('(%3d/%3d) free dof',length(Isup),sum(cg.Nphi));
    statstr = appendstatus({''},str);    
end

% how sparse is each dof
sparseness = Ntot / Npx;
usesparse = mean(sparseness) < corsparseness;

if (cg.precision == 1)
    precision = 'single';
    usesparse = false;
else
    precision = 'double';
end

% memory allocation
U1 = zeros(nroi,mroi,precision);
U2 = zeros(nroi,mroi,precision);
U3 = 0;
U4 = 0;
U5 = 0;
U6 = 0;
if Ndims >= 3
    U3 = zeros(nroi,mroi,precision);
end
if Ndims >= 4
    U4 = zeros(nroi,mroi,precision);
end
if Ndims >= 5
    U5 = zeros(nroi,mroi,precision);
end
if Ndims >= 6
    U6 = zeros(nroi,mroi,precision);
end
if cg.memsave == 0
    % full L matrix definition
    if usegpu
        L = gpuArray.zeros(Npx,Ndof,precision);
    elseif usesparse
        L = spalloc(Npx,Ndof,sum(Ntot));
    else
        L = zeros(Npx,Ndof,precision);
    end
elseif cg.memsave == 1
    % per dimension L matrix definition
    if usegpu
        Li = gpuArray.zeros(Npx,cg.Nphi(1),precision);
        Lj = gpuArray.zeros(Npx,cg.Nphi(1),precision);
    elseif usesparse
        Li = spalloc(Npx,cg.Nphi(1),max(Ntot(Idims{1}))*cg.Nphi(1));
        Lj = spalloc(Npx,cg.Nphi(1),max(Ntot(Idims{1}))*cg.Nphi(1));
    else
        Li = zeros(Npx,cg.Nphi(1),precision);
        Lj = zeros(Npx,cg.Nphi(1),precision);
    end
elseif cg.memsave == 2
    % per dof L matrix definition
    if usegpu
        Li = gpuArray.zeros(Npx,1,precision);
        Lj = gpuArray.zeros(Npx,1,precision);
    elseif usesparse
        Li = spalloc(Npx,1,max(Ntot(Idims{1})));
        Lj = spalloc(Npx,1,max(Ntot(Idims{1})));
    else
        Li = zeros(Npx,1,precision);
        Lj = zeros(Npx,1,precision);
    end
end

% old displacements (to compute the incremental displacement)
Uo1 = zeros(nroi,mroi,precision);
Uo2 = zeros(nroi,mroi,precision);
Uo3 = zeros(nroi,mroi,precision);

% counters
conv = false;
it = 0;
div = 0;

% initialize 
ndu = Inf;
nb = Inf;
ndp = Inf;
dr = -Inf;
convcrit = [ndu ; nb ; ndp ; dr];

convparstr{1,1} = '|du|';
convparstr{2,1} = '|b|';
convparstr{3,1} = '|dp|';
convparstr{4,1} = '|dr|';

str = sprintf('%3s, %10s, %10s, %10s, %10s, %4s, %3s','it',convparstr{1},convparstr{2},convparstr{3},convparstr{4},'div','sp');
statstr = appendstatus(statstr,str);
if headlessmode
    headlessstatus(str);
end

stat = get(S.corstatus,'String');
stat = [sprintf('%2s,%9s,%13s,%4s','it',convparstr{cg.convparam},'mean|r|','div') ; stat];
set(S.corstatus,'String',stat);drawnow
while ~conv
    it = it + 1;
    % Check for stop button
    drawnow
    if get(S.corstop,'Value')
        dic.converged = 4;
        return
    end
    
    if cg.memsave == 0
        if issparse(L)
            sprs = 1;
        else
            sprs = 0;
        end
    else
        if issparse(Li)
            sprs = 1;
        else
            sprs = 0;
        end
    end
    
    % update displacement field
    U1 = reshape(phi12(Iroi,:) * p(Idims{1}),nroi,mroi);
    U2 = reshape(phi12(Iroi,:) * p(Idims{2}),nroi,mroi);
    if Ndims >= 3
        U3 = reshape(phi34(Iroi,:) * p(Idims{3}),nroi,mroi);
    end
    if Ndims >= 4
        U4 = reshape(phi34(Iroi,:) * p(Idims{4}),nroi,mroi);
    end
    if Ndims >= 5
        U5 = reshape(phi34(Iroi,:) * p(Idims{5}),nroi,mroi);
    end
    if Ndims >= 6
        U6 = reshape(phi34(Iroi,:) * p(Idims{6}),nroi,mroi);
    end
    
    % update gtilde
    if usegpu
        gt = gpuArray(interp2(X,Y,gather(g),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0)) ;
        % correct the brightness
        q = gather(U3) + gather(U4).*gt + gather(U5).*gt.^2 + gather(U6).*gt.^3;
    else
        gt = interp2(X,Y,g,X(In,Im)+U1,Y(In,Im)+U2,'spline',0);
        % correct the brightness
        q = U3 + U4.*gt + U5.*gt.^2 + U6.*gt.^3;
    end
    
    if it == 1
        du = [ndu ndu ndu];
    else
        % maximum delta displacements
        du(1) = max(abs(U1(:) - Uo1(:)));
        du(2) = max(abs(U2(:) - Uo2(:)));
        du(3) = max(abs(q(:)  - Uo3(:)));
    end
    
    % update the old displacements
    Uo1 = U1;
    Uo2 = U2;
    Uo3 = q;
    
    % update residual
    r = froi - (gt + q);
    
    % set masked residual to zero
    r(Imask) = 0;
    
    % update the liveview
    if cg.liveview == 2
        lview.it = it;
        lview.r = r;
        lview.U1 = U1(Iq);
        lview.U2 = U2(Iq);
        plotliveview(H,S,lview);
    elseif (cg.liveview == 1) && (it == 1)
        lview.it = it;
        lview.r = r;
        lview.U1 = U1(Iq);
        lview.U2 = U2(Iq);
        plotliveview(H,S,lview);
    end
    
    % update image gradient
    if cg.gradient == 2 %gradfg
        if usegpu
            gradx = 0.5*gpuArray(dfdx + interp2(X,Y,gather(dgdx),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
            grady = 0.5*gpuArray(dfdy + interp2(X,Y,gather(dgdy),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
        else
            gradx = 0.5*(dfdx + interp2(X,Y,dgdx,X(In,Im)+U1,Y(In,Im)+U2,'spline',0));
            grady = 0.5*(dfdy + interp2(X,Y,dgdy,X(In,Im)+U1,Y(In,Im)+U2,'spline',0));
        end
        Lupdate = true;
    elseif cg.gradient == 3 %gradf
        if it == 1
            gradx = dfdx;
            grady = dfdy;
            Lupdate = true;
        else
            Lupdate = false;
        end
    elseif cg.gradient == 4 %gradg
        if usegpu
            gradx = gpuArray(interp2(X,Y,gather(dgdx),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
            grady = gpuArray(interp2(X,Y,gather(dgdy),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
        else
            gradx = interp2(X,Y,dgdx,X(In,Im)+U1,Y(In,Im)+U2,'spline',0);
            grady = interp2(X,Y,dgdy,X(In,Im)+U1,Y(In,Im)+U2,'spline',0);
        end
        Lupdate = true;
    elseif cg.gradient == 1 %auto
        % auto grad threshold
        if (it > 2) && (convcrit(cg.convparam) < corgradcrit*cg.convcrit)
            % close to convergence stop updating the gradient
            Lupdate = false;
        else
            % far from convergence use gradfg
            if usegpu
                gradx = 0.5*gpuArray(dfdx + interp2(X,Y,gather(dgdx),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
                grady = 0.5*gpuArray(dfdy + interp2(X,Y,gather(dgdy),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0));
            else
                gradx = 0.5*(dfdx + interp2(X,Y,dgdx,X(In,Im)+U1,Y(In,Im)+U2,'spline',0));
                grady = 0.5*(dfdy + interp2(X,Y,dgdy,X(In,Im)+U1,Y(In,Im)+U2,'spline',0));
            end
            Lupdate = true;
        end
        
    end
    
    % update L, M and B
    % =================================
    if Lupdate
        if usegpu
            M = gpuArray.zeros(Ndof,Ndof);
            b = gpuArray.zeros(Ndof,1);
        else
            M = zeros(Ndof,Ndof);
            b = zeros(Ndof,1);
        end
        Nsup = zeros(Ndof,1);
        Ntot = zeros(Ndof,1);
    end
    
    % Check for stop button
    drawnow
    if get(S.corstop,'Value')
        dic.converged = 4;
        return
    end
    
    if cg.memsave == 0
        % full L matrix definition (most memory and fastest)

        % update status
        stat = get(S.corstatus,'String');
        stat = ['computing L, M and b (fast)' ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        if Lupdate
            
            % x
            L(:,Idims{1}) = repmat(gradx(Iunmask),1,cg.Nphi(1)) .* phi12(Iroi(Iunmask),:);
            
            % y
            L(:,Idims{2}) = repmat(grady(Iunmask),1,cg.Nphi(2)) .* phi12(Iroi(Iunmask),:);

            % z (or constant relaxation)
            if Ndims >= 3
                L(:,Idims{3}) = phi34(Iroi(Iunmask),:);
            end
            % linear relaxation
            if Ndims >= 4
                L(:,Idims{4}) = repmat(froi(Iunmask),1,cg.Nphi(4)) .* phi34(Iroi(Iunmask),:);
            end
            % quadratic relaxation
            if Ndims >= 5
                L(:,Idims{5}) = repmat(froi(Iunmask).^2,1,cg.Nphi(5)) .* phi34(Iroi(Iunmask),:);
            end
            % cubic relaxation
            if Ndims == 6
                L(:,Idims{6}) = repmat(froi(Iunmask).^3,1,cg.Nphi(6)) .* phi34(Iroi(Iunmask),:);
            end
        end
        
        % update b
        b = L' * r(Iunmask);
        
        % update M
        if Lupdate
            M = L' * L;
        end
        % Check for stop button
        drawnow
        if get(S.corstop,'Value')
            dic.converged = 4;
            return
        end
        
    elseif cg.memsave == 1
        % per dimension L matrix definition (reduced memory by factor Ndims but slower)
        
        % update status
        stat = get(S.corstatus,'String');
        stat = [sprintf('computing L, M and b (%3d/%3d)',0,Ndims) ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        for i = 1:Ndims
            I = Idims{i};
            
            idim = Idofdim(I);
            idof = Idofdof(I);
            
            if idim == 1
                Li(:,1:cg.Nphi(i)) = repmat(gradx(Iunmask),1,cg.Nphi(i)) .* phi12(Iroi(Iunmask),idof);
            elseif idim == 2
                Li(:,1:cg.Nphi(i)) = repmat(grady(Iunmask),1,cg.Nphi(i)) .* phi12(Iroi(Iunmask),idof);
            elseif idim == 3
                Li(:,1:cg.Nphi(i)) = phi34(Iroi(Iunmask),idof);
            elseif idim == 4
                Li(:,1:cg.Nphi(i)) = repmat(froi(Iunmask),1,cg.Nphi(i)) .* phi34(Iroi(Iunmask),idof);
            elseif idim == 5
                Li(:,1:cg.Nphi(i)) = repmat(froi(Iunmask).^2,1,cg.Nphi(i)) .* phi34(Iroi(Iunmask),idof);
            elseif idim == 6
                Li(:,1:cg.Nphi(i)) = repmat(froi(Iunmask).^3,1,cg.Nphi(i)) .* phi34(Iroi(Iunmask),idof);
            end
            
            % update b
            b(I) = Li(:,1:cg.Nphi(i))' * r(Iunmask);
            
            % Check for stop button
            drawnow
            if get(S.corstop,'Value')
                dic.converged = 4;
                return
            end
            
            if Lupdate
                
                for j = 1:Ndims
                    J = Idims{j};
                    
                    jdim = Idofdim(J);
                    jdof = Idofdof(J);
                    
                    if jdim == 1
                        Lj(:,1:cg.Nphi(j)) = repmat(gradx(Iunmask),1,cg.Nphi(j)) .* phi12(Iroi(Iunmask),jdof);
                    elseif jdim == 2
                        Lj(:,1:cg.Nphi(j)) = repmat(grady(Iunmask),1,cg.Nphi(j)) .* phi12(Iroi(Iunmask),jdof);
                    elseif jdim == 3
                        Lj(:,1:cg.Nphi(j)) = phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 4
                        Lj(:,1:cg.Nphi(j)) = repmat(froi(Iunmask),1,cg.Nphi(j)) .* phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 5
                        Lj(:,1:cg.Nphi(j)) = repmat(froi(Iunmask).^2,1,cg.Nphi(j)) .* phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 6
                        Lj(:,1:cg.Nphi(j)) = repmat(froi(Iunmask).^3,1,cg.Nphi(j)) .* phi34(Iroi(Iunmask),jdof);
                    end

                    % update M
                    M(I,J) = Li(:,1:cg.Nphi(i))' * Lj(:,1:cg.Nphi(j));
                end
            end
            
            % update status
            stat{1} = sprintf('computing L, M and b (%3d/%3d)',i,Ndims);
            set(S.corstatus,'String',stat);drawnow
            
            % Check for stop button
            drawnow
            if get(S.corstop,'Value')
                dic.converged = 4;
                return
            end
        end
        
    elseif cg.memsave == 2
        % per dof L matrix definition (reduced memory by factor Ndof but super slow)

        % update status
        stat = get(S.corstatus,'String');
        stat = [sprintf('computing L, M and b (%3d/%3d)',0,Ndof) ; stat];
        set(S.corstatus,'String',stat);drawnow
        
        for i = 1:Ndof
            idim = Idofdim(i);
            idof = Idofdof(i);
            
            if idim == 1
                Li(:,1) = gradx(Iunmask) .* phi12(Iroi(Iunmask),idof);
            elseif idim == 2
                Li(:,1) = grady(Iunmask) .* phi12(Iroi(Iunmask),idof);
            elseif idim == 3
                Li(:,1) = phi34(Iroi(Iunmask),idof);
            elseif idim == 4
                Li(:,1) = froi(Iunmask) .* phi34(Iroi(Iunmask),idof);
            elseif idim == 5
                Li(:,1) = froi(Iunmask).^2 .* phi34(Iroi(Iunmask),idof);
            elseif idim == 6
                Li(:,1) = froi(Iunmask).^3 .* phi34(Iroi(Iunmask),idof);
            end
            
            % update b
            b(i) = Li' * r(Iunmask);
            
            % Check for stop button
            drawnow
            if get(S.corstop,'Value')
                dic.converged = 4;
                return
            end
            
            if Lupdate
                
                for j = 1:Ndof
                    
                    jdim = Idofdim(j);
                    jdof = Idofdof(j);
                    
                    if jdim == 1
                        Lj(:,1) = gradx(Iunmask) .* phi12(Iroi(Iunmask),jdof);
                    elseif jdim == 2
                        Lj(:,1) = grady(Iunmask) .* phi12(Iroi(Iunmask),jdof);
                    elseif jdim == 3
                        Lj(:,1) = phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 4
                        Lj(:,1) = froi(Iunmask) .* phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 5
                        Lj(:,1) = froi(Iunmask).^2 .* phi34(Iroi(Iunmask),jdof);
                    elseif jdim == 6
                        Lj(:,1) = froi(Iunmask).^3 .* phi34(Iroi(Iunmask),jdof);
                    end
                    
                    % update M
                    M(i,j) = Li' * Lj;
                end
            end

            % update status
            stat{1} = sprintf('computing L, M and b (%3d/%3d)',i,Ndof);
            set(S.corstatus,'String',stat);drawnow
            
            % Check for stop button
            drawnow
            if get(S.corstop,'Value')
                dic.converged = 4;
                return
            end
        end
    end
    
    % update status
    stat = stat(2:end);
    set(S.corstatus,'String',stat);drawnow

    % Thikonov regularization
    % =================================
    if Lupdate
        % get the eigenvectors of M
        eigval = eig(M);
    end
    if cg.tikhsteps == 0;
        alpha = 0;
    elseif it <= cg.tikhsteps
        alpha = cg.tikhpar(it)*max(eigval);
    elseif it > cg.tikhsteps
        alpha = cg.tikhpar(end)*max(eigval);
    end
    % create tikhonov matrix and right hand
    Mt = alpha * eye(Ndof);
    bt = alpha * (cg.p - p);
    
    % Solve dp
    % =================================
    if usegpu
        dp = gpuArray.zeros(Ndof,1);
    else
        dp = zeros(Ndof,1);
    end
    % only solve for those dof which have support
    dp(Isup) = (M(Isup,Isup) + Mt(Isup,Isup)) \ (b(Isup) + bt(Isup));
    
    % update p
    p = p + dp;
    
    % store results
    nb = norm(b(Isup) + bt(Isup))/Npx;
    ndp = norm(dp(Isup));
    mr = mean(abs(r(Iunmask))) * 100;
    ndu = norm(du);
    if it == 1
        dr = mr;
    else
        dr = mr - itstore(it-1).mr;
    end
    convcrit = [ndu ; nb ; ndp ; dr];
    
    % check for divergence (dr should always be negative)
    if  (it > 3) && (convcrit(4) > cordivtol)
        div = div + 1;
    else
        div = div - 0.5;
    end
    div = max([div 0]);
    
    % use gather to force all arrays to be CPU (instead of GPU)
    itstore(it).it = gather(it);
    itstore(it).mr = gather(mr);
    itstore(it).convcrit = gather(convcrit);
    itstore(it).p = gather(p);
    itstore(it).div = gather(div);
    itstore(it).tikhonov = gather(alpha);
    
    % status update
    str = sprintf('%3d, %10.3e, %10.3e, %10.3e, %10.3e, %4.1f, %3d',it,convcrit,div,sprs);
    statstr = appendstatus(statstr,str);
    if headlessmode
        headlessstatus(str);
    end
    
    % corstatus update
    stat = get(S.corstatus,'String');
    stat = [sprintf('%2d,%9.2e,%13.6e,%4.1f',it,convcrit(cg.convparam),mr,div) ; stat];
    set(S.corstatus,'String',stat);drawnow
    
    % check for convergence
    if (it > 2) && (abs(convcrit(cg.convparam)) < abs(cg.convcrit))
        converged = 1;
        break
    end
    
    % check for max divergence
    if div >= cg.maxdiv
        converged = 2;
        break
    end
    
    % check for maxit
    if it >= cg.maxit
        converged = 3;
        break
    end

    % Check for continue button
    drawnow
    if get(S.corcontinue,'Value')
        converged = 4;
        set(S.corcontinue,'Value',0)
        drawnow
        break
    end
    
    % Check for stop button
    drawnow
    if get(S.corstop,'Value')
        converged = 5;
        return
    end
end

% debugging
% assignin('base','M',M)
% assignin('base','L',L)


convstr{1,1} = 'converged';
convstr{2,1} = 'diverged';
convstr{3,1} = 'max. it.';
convstr{4,1} = 'continued';
convstr{5,1} = 'stopped';

% choose the best iteration
if strcmp(cg.bestit,'residual')
    bestlist = [itstore.mr];
    [bestlist, I] = sort(bestlist);
    bestit = I(1);
    updatebestit = 1;
elseif strcmp(cg.bestit,'change in disp.')
    bestlist = [itstore.convcrit];
    [bestlist, I] = sort(bestlist(1,:));
    bestit = I(1);
    updatebestit = 1;
elseif strcmp(cg.bestit,'right hand member')
    bestlist = [itstore.convcrit];
    [bestlist, I] = sort(bestlist(2,:));
    bestit = I(1);
    updatebestit = 1;
elseif strcmp(cg.bestit,'update in dof')
    bestlist = [itstore.convcrit];
    [bestlist, I] = sort(bestlist(3,:));
    bestit = I(1);
    updatebestit = 1;
elseif strcmp(cg.bestit,'change in residual')
    bestlist = [itstore.convcrit];
    [bestlist, I] = sort(abs(bestlist(4,:)));
    bestit = I(1);
    updatebestit = 1;
else
    bestit = it;
    updatebestit = 0;
end
    
% update displacement on best it
if updatebestit    
    
    p = itstore(bestit).p;
    
    % update displacement field
    U1 = reshape(phi12(Iroi,:) * p(Idims{1}),nroi,mroi);
    U2 = reshape(phi12(Iroi,:) * p(Idims{2}),nroi,mroi);
    if Ndims >= 3
        U3 = reshape(phi34(Iroi,:) * p(Idims{3}),nroi,mroi);
    end
    if Ndims >= 4
        U4 = reshape(phi34(Iroi,:) * p(Idims{4}),nroi,mroi);
    end
    if Ndims >= 5
        U5 = reshape(phi34(Iroi,:) * p(Idims{5}),nroi,mroi);
    end
    if Ndims == 6
        U6 = reshape(phi34(Iroi,:) * p(Idims{6}),nroi,mroi);
    end
    
    % update gtilde
    if usegpu
        gt = gpuArray(interp2(X,Y,gather(g),X(In,Im)+gather(U1),Y(In,Im)+gather(U2),'spline',0)) ;
        % correct the brightness
        q = gather(U3) + gather(U4).*gt + gather(U5).*gt.^2 + gather(U6).*gt.^3;
    else
        gt = interp2(X,Y,g,X(In,Im)+U1,Y(In,Im)+U2,'spline',0);
        % correct the brightness
        q = U3 + U4.*gt + U5.*gt.^2 + U6.*gt.^3;
    end

    % update residual
    r = froi - (gt + q);
    r(Imask) = 0;
end

% test memory availability
if usegpu
    D.gpu = gpuDevice;
    free = D.gpu.FreeMemory;    
else
    free = memfree;
end

% status
str = sprintf('%s, using iteration %d (free memory: %g Mb)',convstr{converged},bestit,round(free*1e-6));
statstr = appendstatus(statstr,str);
if headlessmode
    headlessstatus(str);
end


stat = get(S.corstatus,'String');
stat = [sprintf('%s, using iteration %d',convstr{converged},bestit) ; stat];
set(S.corstatus,'String',stat);drawnow

% update final liveview
if cg.liveview >= 1
    lview.it = it;
    lview.r = r;
    lview.U1 = U1(Iq);
    lview.U2 = U2(Iq);
    plotliveview(H,S,lview);
end

% info table
Nit = it;
itable = zeros(Nit,9);
for k = 1:Nit
    itable(k,1) = cg.inc;
    itable(k,2) = cg.icg;
    itable(k,3) = k;
    itable(k,4) = Ndof;
    itable(k,5) = Npx;
    itable(k,6) = itstore(k).mr;
    itable(k,7:10) = itstore(k).convcrit;
    itable(k,11) = itstore(k).tikhonov;
    itable(k,12) = itstore(k).div;
    if k == bestit;
        itable(k,13) = converged;
    else
        itable(k,13) = 0;
    end
end

if cg.dimensions == 3 % 3D
    % define the z displacment as positive to the top
    U3 = -U3;
end

dic.itstore = itstore;
dic.itable = itable;
dic.p = p;
dic.U1 = U1;
dic.U2 = U2;
if Ndims >= 3
    dic.U3 = U3;
end
if Ndims >= 4
    dic.U4 = U4;
end
if Ndims >= 5
    dic.U5 = U5;
end
if Ndims == 6
    dic.U6 = U6;
end
dic.froi = froi;
dic.gt = gt;
dic.xroi = xroi;
dic.yroi = yroi;
dic.r = r;
dic.q = q;
dic.In = In;
dic.Im = Im;
dic.Imask = Imask;
dic.Iunmask = Iunmask;
dic.statstr = statstr;
dic.converged = converged;




function [] = plotliveview(H,S,lview)
D = guidata(H);

% assignin('base','lview',lview);

% plot image
set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes7)
colormap(D.gui.color.cmap)
title(sprintf('residual [%%], inc %d, cg %d, it %d',lview.inc,lview.icg,lview.it))
hi = findobj(H,'Tag','cor_image');
if isempty(hi)
    set(gca,'NextPlot','replacechildren')
    imagesc(lview.xroi,lview.yroi,lview.r*100,'Parent',S.axes7,'cor_image');
    set(gca,'color','k')
    xlabel('x [px]')
    ylabel('y [px]')
    set(gca,'ydir','reverse')
    daspect([1 1 1 ])
    hc = colorbar;
    set(hc,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'FontSize',D.gui.fontsize(3))
else
    set(hi,'CData',lview.r*100,'XData',lview.xroi([1,end]),'YData',lview.yroi([1,end]));
end
set(gca,'xlim',lview.xlim)
set(gca,'ylim',lview.ylim)

clim = mean(lview.r(:)) + [-3 3]*std(lview.r(:));
clim = clim*100;
if clim(1) == clim(2)
    clim = clim(1) + abs(clim(1))*[-0.01 0.01];
end

set(gca,'NextPlot','add');
hq = findobj(H,'Tag','cor_quiver');
if isempty(hq)
    hq = quiver(lview.X,lview.Y,lview.U1,lview.U2,0,'Tag','cor_quiver');
    set(hq,'color',D.gui.color.fg);
else
    set(hq,'XData',lview.X,'YData',lview.Y,'UData',lview.U1,'VData',lview.U2);
end

ht = findobj(H,'Tag','cor_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',lview.inc,lview.Ninc),'units','normalized','Tag','cor_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(3))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',lview.inc,lview.Ninc));
end

caxis(gather(clim));
drawnow



