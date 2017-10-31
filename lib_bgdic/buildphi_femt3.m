function phi = buildphi_femt3(x,y,phi_crds,phi_conn,H,varargin)
%BUILDPHI_FEMT3 creates the basis function matrix phi, which has one row for
%    each pixel and one column for each basis function.
%
%    phi = buildphi_femt3(x,y,phi_crds,phi_conn)
%
%    INPUTS:
%    x,y      : the image coordinate vectors of size m and n respectivily
%    phi_crds : a (N x 2) matrix with node coordinates
%               the number of basis functions will be equal to N
%    phi_conn : a (M x 3) matrix with the connectivity of each element
%             
%
%    OPTIONAL INPUTS:
%    phi = buildphi_femt3(x,y,phi_crds,phi_conn,'single')
%    By default BUILDPHI_FEMT3 uses double precision floating point values
%    for the the phi matrix, however when running into memory issues it is
%    possible to use singe precision, nevertheless this will influence the
%    DIC accuracy greatly and is not recommended.
%
%    OUTPUT:
%    Only one version of the basis function matrix is stored, i.e. phi.a.
%    This is to same memory because phi.x, phi.y, phi.z would be exactly
%    equal.
%
%    EXAMPLE:
%    x = linspace(-5,5,50);
%    y = linspace(-4,4,40);
%    phi_crds(1,:) = [-4 -3];
%    phi_crds(2,:) = [ 4 -3];
%    phi_crds(3,:) = [ 4  3];
%    phi_crds(4,:) = [-4  3];
%    phi_crds(5,:) = [ 0  0];
%    phi_conn(1,:) = [ 1 2 5];
%    phi_conn(2,:) = [ 2 3 5];
%    phi_conn(3,:) = [ 3 4 5];
%    phi_conn(4,:) = [ 4 1 5];
%    phi = buildphi_femt3(x,y,phi_crds,phi_conn);
%    phiplot(x,y,phi);
%
%    See also JNDIC, PHIPLOT, BUILDPHI_POLY, BUILDPHI_LEGENDRE, 
%             BUILDPHI_CHEBYSHEV.
%
%copyright: Jan Neggers, 2013

precision = 'double';
if nargin == 6
    if strcmpi(varargin{1},'single')
        precision = 'single';
    elseif strcmpi(varargin{1},'double')
        precision = 'double';
    else
        error('buildphi_femt3: unknown precision options')
    end
end

nc = phi_crds;
ec = phi_conn;

[X, Y] = meshgrid(x,y);
[n, m] = size(X);
N   = size(nc,1);
Nel =  size(ec,1);

% initiate shape functions
phi.a = zeros(n*m,N,precision);

% store image size
phi.n = n;
phi.m = m;

% element node selectors
p = 1;
q = 2;
r = 3;
% loop over each element
for k = 1:Nel
    % find the pixels which are inside the element
    In = find(inpolygon(X,Y,nc(ec(k,:),1),nc(ec(k,:),2)));
    
    for kk = 1:3
        % select the node coordinates
        xp = nc( ec(k,p) ,1);
        yp = nc( ec(k,p) ,2);
        xq = nc( ec(k,q) ,1);
        yq = nc( ec(k,q) ,2);
        xr = nc( ec(k,r) ,1);
        yr = nc( ec(k,r) ,2);
        
        % calculate the plane
        a = (yr-yq)/(xr*yq-xq*yr-xr*yp+xp*yr+xq*yp-xp*yq);
        b = -(xr-xq)/(xr*yq-xq*yr-xr*yp+xp*yr+xq*yp-xp*yq);
        c = (xr*yq-xq*yr)/(xr*yq-xq*yr-xr*yp+xp*yr+xq*yp-xp*yq);
        N = a.*X+b.*Y+c;
        
        % get the current node number
        node = ec(k,p);

        % occasionally inpolygon counts edge pixels twice
        I = find(phi.a(:,node)~=0);
        
        % find all pixels that do not have a shapefun value
        In = setdiff(In,I);
        
        % store in the phi.x matrix
        phi.a(In,node) = phi.a(In,node) + N(In);

        % rotate the element node selectors
        [p, q, r] = deal(q,r,p);
    end
    if mod(k,floor(Nel/10)) == 1
        bcwaitbar(H,k/Nel,sprintf('building T3 basis (%d/%d)',k,Nel));
    end
end
bcwaitbar(H);

% partition of unity test
% phi.x(:,5) = sum(phi.x,2);
