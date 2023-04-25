function phi = buildphi_bspline(x,y,xknot,yknot,p,roi,H,varargin)
%BUILDPHI_POLY creates the basis function matrix phi, which has one row for
%    each pixel and one column for each basis function.
%
%    buildphi_bspline(x,y,xknots,yknots,p,roi)
%
%
%copyright: Jan Neggers, 2013

precision = 'double';
if nargin == 8
    if strcmpi(varargin{1},'single')
        precision = 'single';
    elseif strcmpi(varargin{1},'double')
        precision = 'double';
    else
        error('buildphi_bspline: unknown precision options')
    end
end

if false
    p = mesh.p;
    xknot = mesh.xknot;
    yknot = mesh.yknot;
    load('mesh.mat','mesh');
    roi = [51.2000  460.8000   38.4000  345.6000];
    n = 384;
    m = 512;
    x = 1:m;
    y = 1:n;
end

n = length(y);
m = length(x);

if length(roi) == 1
    roi = [min(x)+roi max(x)-roi min(y)+roi max(y)-roi];
end

Im = (x < roi(1) | x > roi(2));
In = (y < roi(3) | y > roi(4));

% build phi
% ======================
Px = bsplines1D(x,xknot,p);
Py = bsplines1D(y,yknot,p);
% zero outside roi
Px(Im,:) = 0;
Py(In,:) = 0;

% number of dof in each direction
Nx = size(Px,2);
Ny = size(Py,2);

% initiate phi matrix
Ndof = Nx*Ny;
Npx = n*m;
phi = zeros(Npx,Ndof);

% dyadic multiplication
k = 0;
for kx = 1:Nx
    for ky = 1:Ny
        k = k + 1;
        P = Py(:,ky) * Px(:,kx)';
        phi(:,k) = P(:);
    end
    bcwaitbar(H,kx/Nx,sprintf('building B-Spline basis (%d/%d)',k,Nx*Ny));
end
bcwaitbar(H);

if strcmpi(precision,'single')
    phi = single(phi);
end

