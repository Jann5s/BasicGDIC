function phi = buildphi_poly(x,y,phi_list,roi,varargin)
%BUILDPHI_POLY creates the basis function matrix phi, which has one row for
%    each pixel and one column for each basis function.
%
%    buildphi_poly(x,y,phi_list,roi)
%
%    INPUTS:
%    n,m      : the image dimentions of the ROI, the shape functions will
%               be cacluated on a domain -1 <= x,y < 1 with n steps in y
%               direction and m steps in x direction
%    phi_list : a list of basis functions of size (N x 3) where N is the
%               number of basis functions and each row contains the three
%               parameters, [a b c], where [a b] are the polynomial orders:
%               P = (x^a)*(y^b), and [c] is the direction in which the basis
%               function is applied, 1=x, 2=y, and 3=z.
%
%    OPTIONAL INPUTS:
%    phi = buildphi_poly(x,y,phi_list,roi,'single')
%    By default BUILDPHI_POLY uses double precision floating point values
%    for the the phi matrix, however when running into memory issues it is
%    possible to use singe precision, nevertheless this will influence the
%    DIC accuracy greatly and is not recommended.
%
%    OUTPUT:
%    The basis function matrix is split into a parts according to the
%    application direction, i.e. phi.x, phi.y, phi.z. In other words, this
%    function returns a structure with one matrix per dimension.
%
%    EXAMPLE:
%    roi = [-5 5 -4 4];
%    [n m] = deal(120,150);
%    phi_list(1,:) = [0 0 1];  % x-translation
%    phi_list(2,:) = [0 0 2];  % y-translation
%    phi_list(3,:) = [0 0 3];  % z-translation
%    phi_list(4,:) = [0 1 1];  % strain-xx
%    phi_list(5,:) = [0 1 2];  % shear-xy
%    phi_list(6,:) = [0 1 3];  % tilt around y-axis
%    phi_list(7,:) = [1 0 1];  % shear-yx
%    phi_list(8,:) = [1 0 2];  % strain-yy
%    phi_list(9,:) = [1 0 3];  % tilt around x-axis
%    phi = buildphi_poly(x,y,phi_list,roi);
%    phiplot(x,y,phi);
%
%    See also JNDIC, PHIPLOT, BUILDPHI_LEGENDRE, BUILDPHI_CHEBYSHEV, 
%             BUILDPHI_FEMT3.
%
%copyright: Jan Neggers, 2013

precision = 'double';
if nargin == 5
    if strcmpi(varargin{1},'single')
        precision = 'single';
    elseif strcmpi(varargin{1},'double')
        precision = 'double';
    else
        error('buildphi_poly: unknown precision options')
    end
end

% Debugging
% ==================
if ~true
    clear all
    precision = 'double';
    roi = [-5 5 -4 4];
    [n m] = deal(120,150);
    x = linspace(roi(1)-2,roi(2)+2,m+2); x = x(2:end-1);
    y = linspace(roi(3)-2,roi(4)+2,n+2); y = y(2:end-1);
    phi_list(1,:) = [0 0 1];  % x-translation
    phi_list(2,:) = [0 0 2];  % y-translation
    phi_list(3,:) = [0 0 3];  % z-translation
    phi_list(4,:) = [0 1 1];  % strain-xx
    phi_list(5,:) = [0 1 2];  % shear-xy
    phi_list(6,:) = [0 1 3];  % tilt around y-axis
    phi_list(7,:) = [1 0 1];  % shear-yx
    phi_list(8,:) = [1 0 2];  % strain-yy
    phi_list(9,:) = [1 0 3];  % tilt around x-axis
    phi_list(10,:) = [0 2 1];
    phi_list(11,:) = [0 2 2];
    phi_list(12,:) = [0 2 3];
    phi_list(13,:) = [1 1 1];
    phi_list(14,:) = [1 1 2];
    phi_list(15,:) = [1 1 3];
    phi_list(16,:) = [2 0 1];
    phi_list(17,:) = [2 0 2];
    phi_list(18,:) = [2 0 3];
    
    % Remove DOF in y-direction
    phi_list = phi_list(phi_list(:,3)~=2,:);
    
%     phi_list = [1 2 3];
end

if length(roi) == 1
    roi = [min(x)+roi max(x)-roi min(y)+roi max(y)-roi];
end

Im = find(x < roi(1) | x > roi(2));
In = find(y < roi(3) | y > roi(4));

% create a space
x = 2*(x - 0.5*(roi(1)+roi(2))) / (roi(2) - roi(1));
y = 2*(y - 0.5*(roi(3)+roi(4))) / (roi(4) - roi(3));

m = length(x);
n = length(y);

% split list into directions
listx = phi_list(phi_list(:,3)==1,:);
listy = phi_list(phi_list(:,3)==2,:);
listz = phi_list(phi_list(:,3)==3,:);

% number of DOF per direction
Nx = size(listx,1);
Ny = size(listy,1);
Nz = size(listz,1);

% initiate matrices
phi.x = zeros(n*m,Nx,precision);
phi.y = zeros(n*m,Ny,precision);
phi.z = zeros(n*m,Nz,precision);

% store image size
phi.n = n;
phi.m = m;

% build phi.x
% ======================
for k = 1:Nx
    P = y.^listx(k,2) * x.^listx(k,1);
%     P(In,:) = 0;
%     P(:,Im) = 0;
    phi.x(:,k) = P(:);
end

% build phi.y
% ======================
for k = 1:Ny
    P = y.^listy(k,2) * x.^listy(k,1);
%     P(In,:) = 0;
%     P(:,Im) = 0;
    phi.y(:,k) = P(:);
end

% build phi.z
% ======================
for k = 1:Nz
    P = y.^listz(k,2) * x.^listz(k,1);
%     P(In,:) = 0;
%     P(:,Im) = 0;
    phi.z(:,k) = P(:);
end
