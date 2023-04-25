function phi = buildphi_chebyshev(x,y,phi_list,roi,varargin)
%BUILDPHI_CHEBYSHEV creates the basis function matrix phi, which has one 
%    row for each pixel and one column for each basis function. The created
%    2D basis functions are based on the Chebyshev polynomials of the first
%    kind by Pafnuty Chebyshev (alternatively transliterated as Chebychev,
%    Chebysheff, Chebyshov, Tchebychev or Tchebycheff, Tschebyschev or 
%    Tschebyscheff).
%
%    phi = buildphi_chebyshev(x,y,phi_list,roi)
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
%    phi = buildphi_chebyshev(x,y,phi_list,roi,'single')
%    By default BUILDPHI_CHEBYSHEV uses double precision floating point
%    values for the the phi matrix, however when running into memory issues
%    it is possible to use singe precision, nevertheless this will
%    influence the DIC accuracy greatly and is not recommended.
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
%    phi = buildphi_chebyshev(x,y,phi_list,roi);
%    phiplot(x,y,phi);
%
%    See also JNDIC, PHIPLOT, BUILDPHI_POLY, BUILDPHI_LEGENDRE, 
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
        error('buildphi_chebyshev: unknown precision options')
    end
end

if length(roi) == 1
    roi = [min(x)+roi max(x)-roi min(y)+roi max(y)-roi];
end

% create a space
x = 2*(x - 0.5*(roi(1)+roi(2))) / (roi(2) - roi(1));
y = 2*(y - 0.5*(roi(3)+roi(4))) / (roi(4) - roi(3));

m = length(x);
n = length(y);

% number of shape functions
N = size(phi_list,1);

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
    a = listx(k,1);
    b = listx(k,2);
    
    Tx = chebyshev1D(x,a);
    Ty = chebyshev1D(y,b);
%     Tx(Im) = 0;
%     Ty(In) = 0;
    
    P = Ty' * Tx;
    phi.x(:,k) = P(:);
end

% build phi.y
% ======================
for k = 1:Ny
    a = listy(k,1);
    b = listy(k,2);
    
    Tx = chebyshev1D(x,a);
    Ty = chebyshev1D(y,b);
%     Tx(Im) = 0;
%     Ty(In) = 0;
    
    P = Ty' * Tx;
    phi.y(:,k) = P(:);
end

% build phi.z
% ======================
for k = 1:Nz
    a = listz(k,1);
    b = listz(k,2);
    
    Tx = chebyshev1D(x,a);
    Ty = chebyshev1D(y,b);
%     Tx(Im) = 0;
%     Ty(In) = 0;
    
    P = Ty' * Tx;
    phi.z(:,k) = P(:);
end

function T = chebyshev1D(x,p)
% this function builds the 1D polynomial of order p for n points
    
n = length(x);

% initiate the first two polynomials
Tn = ones(1,n);
Tk = x;

% loop until the required order
for a = 0:p
    if a == 0
        T = Tn;
    elseif a == 1
        T = Tk;
    else
        % calculate the new polynomial
        T = 2*x.*Tk - Tn;
        
        % update the previous two
        Tn = Tk;
        Tk = T;
    end
end
