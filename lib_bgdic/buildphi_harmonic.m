function phi = buildphi_harmonic(x,y,p,roi,H,varargin)


precision = 'double';
if nargin == 6
    if strcmpi(varargin{1},'single')
        precision = 'single';
    elseif strcmpi(varargin{1},'double')
        precision = 'double';
    else
        error('buildphi_bspline: unknown precision options')
    end
end

% for some (yet unknown) reason, the harmonic functions work better when
% this factor is 1.5 (such that a the waves fit 1.5*k times the ROI.
% fudgefactor = 1.5;
fudgefactor = 1;

% create a space
x = fudgefactor*2*(x - 0.5*(roi(1)+roi(2))) / (roi(2) - roi(1));
y = fudgefactor*2*(y - 0.5*(roi(3)+roi(4))) / (roi(4) - roi(3));

% x = 1.5*(x - roi(1)) / (roi(2) - roi(1));
% y = 1.5*(y - roi(3)) / (roi(4) - roi(3));

m = length(x);
n = length(y);

% build phi
% ======================
Px = harmonic1D(x,p);
Py = harmonic1D(y,p);
% zero outside roi
% Px(Im,:) = 0;
% Py(In,:) = 0;

% number of dof in each direction
Nx = size(Px,2);
Ny = size(Py,2);

% initiate phi matrix
Ndof = Nx*Ny;
Npx = n*m;
phi = zeros(Npx,Ndof,precision);

% dyadic multiplication
k = 0;
for kx = 1:Nx
    for ky = 1:Ny
        k = k + 1;
        P = Py(:,ky) * Px(:,kx)';
        phi(:,k) = P(:);
    end
    bcwaitbar(H,kx/Nx,sprintf('building harmonic basis (%d/%d)',k,Nx*Ny));
end

bcwaitbar(H);
