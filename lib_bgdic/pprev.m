function p = pprev(H,cor,cg,phi12,phi34)
% Convert the DOF of the previous increment or coarse grain step to the
% current. This requires the basis functions of both the previous and the
% current step and the DOF of the previous step

D = guidata(H);

% the old region of interest
xroi = cor.xroi;
yroi = cor.yroi;

% get the mask defined on the old region of interest
Iunmask = cor.Iunmask;
if isempty(Iunmask)
    warning('empty Iunmask: all pixels are masked, possible bug detected');
end

% number of pixels in previous step
npx = length(xroi)*length(yroi);

% fit x and y direction
% ==================================

% create temporary phi matrix
phi = zeros(npx,cg.Nphi(1));

% interpolate new phi matrix on old coordinates
for i = 1:cg.Nphi(1)
    P = reshape(phi12(:,i),cg.n,cg.m);
    P = interp2(cg.x,cg.y,P,xroi,yroi','linear',NaN);
    phi(:,i) = P(:);
end

% check for unsupported dof
Ntot = sum(phi~=0);
Nsup = sum(phi(Iunmask,:)~=0);

% indices of the well supported shape functions
cordofsupport = D.gui.ini.cordofsupport.value;
Isup = find(Nsup./Ntot > cordofsupport);
Iunsup = find(Nsup./Ntot <= cordofsupport);

p1 = zeros(cg.Nphi(1),1);
p2 = zeros(cg.Nphi(1),1);
% fit new shapefunctions on old displacement field
p1(Isup) = phi(Iunmask,Isup) \ cor.U1(Iunmask);
p2(Isup) = phi(Iunmask,Isup) \ cor.U2(Iunmask);
% combine DOF
p = [p1 ; p2];

% fit 3rd and 4th direction (z or brightness and contrast)
% ==================================
if cg.Ndims >= 3
    % create temporary phi matrix
    phi = zeros(npx,cg.Nphi(3));
    
    % interpolate new phi matrix on old coordinates
    for i = 1:cg.Nphi(3)
        P = reshape(phi34(:,i),cg.n,cg.m);
        P = interp2(cg.x,cg.y,P,xroi,yroi','linear',NaN);
        phi(:,i) = P(:);
    end
    
    % check for unsupported dof
    Ntot = sum(phi~=0);
    Nsup = sum(phi(Iunmask,:)~=0);
    
    % indices of the well supported shape functions
    cordofsupport = D.gui.ini.cordofsupport.value;
    Isup = find(Nsup./Ntot > cordofsupport);
    Iunsup = find(Nsup./Ntot <= cordofsupport);
    
    % fit new shapefunctions on old displacement field
    if isfield(cor,'U3') && ~isempty(cor.U3)
        if cg.dimensions == 3 % 3D
            % revert the correction of the sign for U3
            U3 = -cor.U3;
        else
            U3 = cor.U3;
        end
        p3 = zeros(cg.Nphi(3),1);
        p3(Isup) = phi(Iunmask,Isup) \ U3(Iunmask);
    else
        p3 = zeros(cg.Nphi(3),1);
    end
    p = [p1 ; p2 ; p3];
end
if cg.Ndims >= 4
    
    % fit new shapefunctions on old displacement field
    if isfield(cor,'U4') && ~isempty(cor.U4)
        p4 = zeros(cg.Nphi(4),1);
        p4(Isup) = phi(Iunmask,Isup) \ cor.U4(Iunmask);
    else
        p4 = zeros(cg.Nphi(4),1);
    end
    p = [p1 ; p2 ; p3 ; p4];
end
if cg.Ndims >= 5
    
    % fit new shapefunctions on old displacement field
    if isfield(cor,'U5') && ~isempty(cor.U5)
        p5 = zeros(cg.Nphi(5),1);
        p5(Isup) = phi(Iunmask,Isup) \ cor.U5(Iunmask);
    else
        p5 = zeros(cg.Nphi(5),1);
    end
    p = [p1 ; p2 ; p3 ; p4 ; p5];
end
if cg.Ndims >= 6
    
    % fit new shapefunctions on old displacement field
    if isfield(cor,'U6') && ~isempty(cor.U6)
        p6 = zeros(cg.Nphi(6),1);
        p6(Isup) = phi(Iunmask,Isup) \ cor.U6(Iunmask);
    else
        p6 = zeros(cg.Nphi(6),1);
    end
    p = [p1 ; p2 ; p3 ; p4 ; p5 ; p6];
end