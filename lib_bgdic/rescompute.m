% ==================================================
function [] = rescompute(varargin)
% compute result fields
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% if D.usegpu
%     D.gpu = gpuDevice;
%     guidata(H,D);
% end

if ~isfield(D,'files')
    return;
end
if ~isfield(D,'cor')
    return;
end

% get the popupmenu values
str = get(S.resstraindef,'String');
straindef = str{get(S.resstraindef,'Value')};

% quiver positions
Ncg = 4;
Ninc = length(D.cor);
for inc = 1:Ninc
    if D.cor(inc).done ~= Ncg;
        continue
    end
    
    dx = mean(diff(D.cor(inc).xroi));
    dy = mean(diff(D.cor(inc).yroi));
    
    % get data to cpu (if using gpu)
    U1 = gather(D.cor(inc).U1);
    U2 = gather(D.cor(inc).U2);
    
    % use the mask to set some NaN's
    Imask = D.cor(inc).Imask;
    U1(Imask) = NaN;
    U2(Imask) = NaN;
        
    % deformation gradient tensor
    [Fxx, Fxy] = gradient(U1,dx,dy);
    [Fyx, Fyy] = gradient(U2,dx,dy);
    Fxx = Fxx + 1;
    Fyy = Fyy + 1;
    [n, m] = size(Fxx);
    
    ZERO = zeros(n,m,'int8');
    
    % Cauchy Green deformation tensor
    Cxx = Fxx .* Fxx + Fyx .* Fyx;
    Cxy = Fxx .* Fxy + Fyx .* Fyy;
    Cyx = Fyy .* Fyx + Fxy .* Fxx;
    Cyy = Fyy .* Fyy + Fxy .* Fxy;
    
    % Finger deformation tensor
    Bxx = Fxx .* Fxx + Fxy .* Fxy;
    Bxy = Fxx .* Fyx + Fxy .* Fyy;
    Byx = Fyx .* Fxx + Fyy .* Fxy;
    Byy = Fyx .* Fyx + Fyy .* Fyy;
    
    if strcmp(straindef,'none')
        % Do not compute strain
        Exx = ZERO;
        Eyy = ZERO;
        Exy = ZERO;
        Eyx = ZERO;
    elseif strcmp(straindef,'membrane strain')
        % Membrane Strain
        
        % assuming flat initial situation
        [X, Y] = meshgrid(D.cor(inc).xroi,D.cor(inc).yroi);
        [n, m] = size(X);
        Z = zeros(n,m);
       
        if ~isfield(D.cor(inc),'U3') || isempty(D.cor(inc).U3)
            U3 = Z;
        else
            U3 = gather(D.cor(inc).U3);
        end
        U3(Imask) = NaN;
        
        [dXx, dXy] = gradient(X+U1,dx,dy);
        [dYx, dYy] = gradient(Y+U2,dx,dy);
        [dZx, dZy] = gradient(Z+U3,dx,dy);
        
        % this strain definition works for stretched membranes, but is far from
        % universal, test, check, and verify, before using.
        Exx  = hypot(dXx,dZx) - 1;
        Eyy  = hypot(dYy,dZy) - 1;
        Exy = hypot(dXy,dZx);
        Eyx = hypot(dYx,dZy);
        
    elseif strcmp(straindef,'small strain')
        Exx = Fxx - 1;
        Eyy = Fyy - 1;
        Exy = 0.5*(Fxy + Fyx);
        Eyx = 0.5*(Fxy + Fyx);
    elseif strcmp(straindef,'logarithmic strain')
        C = {Cxx,Cxy;Cyx,Cyy};
        [v, d] = eig2d(C);
        Exx = log(sqrt(d{1})) .* v{1,1} .* v{1,1} ...
            + log(sqrt(d{2})) .* v{2,1} .* v{2,1} ;
        Exy = log(sqrt(d{1})) .* v{1,1} .* v{1,2} ...
            + log(sqrt(d{2})) .* v{2,1} .* v{2,2} ;
        Eyx = log(sqrt(d{1})) .* v{1,2} .* v{1,1} ...
            + log(sqrt(d{2})) .* v{2,2} .* v{2,1} ;
        Eyy = log(sqrt(d{1})) .* v{1,2} .* v{1,2} ...
            + log(sqrt(d{2})) .* v{2,2} .* v{2,2} ;
        
    elseif strcmp(straindef,'Green-Lagrange strain')
        % Green Lagrange Strain tensor
        Exx = 0.5 * (Cxx - 1);
        Exy = 0.5 * (Cxy);
        Eyx = 0.5 * (Cyx);
        Eyy = 0.5 * (Cyy - 1);
    elseif strcmp(straindef,'Euler-Almansi strain')
        B = {Bxx,Bxy;Byx,Byy};
        [v, d] = eig2d(B);
        Exx = (1./d{1}) .* v{1,1} .* v{1,1} ...
            + (1./d{2}) .* v{2,1} .* v{2,1} ;
        Exy = (1./d{1}) .* v{1,1} .* v{1,2} ...
            + (1./d{2}) .* v{2,1} .* v{2,2} ;
        Eyx = (1./d{1}) .* v{1,2} .* v{1,1} ...
            + (1./d{2}) .* v{2,2} .* v{2,1} ;
        Eyy = (1./d{1}) .* v{1,2} .* v{1,2} ...
            + (1./d{2}) .* v{2,2} .* v{2,2} ;
        
        Exx = 0.5*(1-Exx);
        Exy = 0.5*(0-Exy);
        Eyx = 0.5*(0-Eyx);
        Eyy = 0.5*(1-Eyy);
    end
    
    if ~strcmp(straindef,'none')
        % major and minor strain
        E = {Exx,Exy;Eyx,Eyy};
        [v, d] = eig2d(E);
        Q11 = v{1,1};
        Q12 = v{1,2};
        Q21 = v{2,1};
        Q22 = v{2,2};
        Emin = d{1,1};
        Emaj = d{2,1};
        Eeq = sqrt(Emin.^2 + Emaj.^2);
    else
        Emaj = ZERO;
        Emin = ZERO;
        Eeq = ZERO;
        Q11 = ZERO;
        Q12 = ZERO;
        Q21 = ZERO;
        Q22 = ZERO;
    end
    
    % use the mask to set some NaN's
    Imask = D.cor(inc).Imask;
    Exx(Imask) = NaN;
    Eyy(Imask) = NaN;
    Exy(Imask) = NaN;
    Eyx(Imask) = NaN;
    Emaj(Imask) = NaN;
    Emin(Imask) = NaN;
    Eeq(Imask) = NaN;
    
    res(inc).Exx = Exx;
    res(inc).Eyy = Eyy;
    res(inc).Exy = Exy;
    res(inc).Eyx = Eyx;
    res(inc).Emaj = Emaj;
    res(inc).Emin = Emin;
    res(inc).Eeq = Eeq;
    res(inc).Q11 = Q11;
    res(inc).Q12 = Q12;
    res(inc).Q21 = Q21;
    res(inc).Q22 = Q22;
    
    % update status
    stat = sprintf('[8] Results recomputed for increment %d',inc);
    D.gui.stat = appendstatus(D.gui.stat,stat);

    bcwaitbar(H,inc/Ninc,sprintf('computing strain (%d/%d)',inc,Ninc));

end

D.res = res;

% update the application data
guidata(H,D);
bcwaitbar(H);
