function phi = buildphi_femt6(x,y,node,conn,H,varargin)
% phi = buildphi_femt6(x,y,node,conn)
%
% Function for creating computing the basis function values for each pixel
% for each degree of freedom (node), for a mesh with 6 noded triangular
% elements (aka T6).
%
% x (1,m) and y (1,n) are the pixel positions
% node is a matrix with nodal coordinates (Ndof,2)
% conn is a matrix with element connectivity (Nel,6)
%
% the node numbering convention used is
%           3
%         /   \
%        6      5
%      /         \
%     1 --- 4 ---- 2

% use only the interpolation method, or also perform an iterative scheme
% per pixel to optimze the revese mapping.
itsolve = false;

precision = 'double';
if nargin == 6
    if strcmpi(varargin{1},'single')
        precision = 'single';
    elseif strcmpi(varargin{1},'double')
        precision = 'double';
    else
        error('buildphi_femt6: unknown precision options')
    end
end

[X, Y] = meshgrid(x,y);
[n, m] = size(X);

Nnode = size(node,1);
Nel = size(conn,1);

phi = zeros(n*m,Nnode,precision);

convcrit = 1e-6;
maxit = 20;
boundrange = 0.1; % evaluation range outside node polygon (relative to largest nodal distance within element)
tol = 1e-6; % pixel position tolerance

% lambda: area coordinates (barycentric) (note: l1 + l2 + l3 = 1)
% zeta:   isoparametric coordinates (carthesian)

% lambda  as a function of zeta
L = @(zeta) [ (1/3) - zeta(1,:) - (1/sqrt(3))*zeta(2,:);
              (1/3) + zeta(1,:) - (1/sqrt(3))*zeta(2,:);
              (1/3) + (2/sqrt(3))*zeta(2,:)];

% zeta as a function of lambda
Z = @(lambda) [ lambda(2,:) - 0.5 + 0.5*lambda(3,:) ; % alt: 0.5 - lambda(1,:) - 0.5*lambda(3,:)
                0.5*sqrt(3)*(lambda(3,:) - (1/3))   ];

% T6 shape functions
N = @(lambda) [lambda(1,:).*(2*lambda(1,:) - 1);
               lambda(2,:).*(2*lambda(2,:) - 1);
               lambda(3,:).*(2*lambda(3,:) - 1);
               4*lambda(1,:).*lambda(2,:);
               4*lambda(2,:).*lambda(3,:);
               4*lambda(3,:).*lambda(1,:)];

% T6 Derivatives
dN = @(zeta) [ 4*zeta(1,:) + (4*sqrt(3)*zeta(2,:))./3 - 1/3                ,... %dN1/dzeta(1,:)
              (4*zeta(2,:))./3 + (4*sqrt(3)*zeta(1,:))./3 - sqrt(3)/9      ;    %dN1/dzeta(2,:)
               4*zeta(1,:) - (4*sqrt(3)*zeta(2,:))./3 + 1/3                ,... %dN2/dzeta(1,:)
              (4*zeta(2,:))./3 - (4*sqrt(3)*zeta(1,:))./3 - sqrt(3)/9      ;    %dN2/dzeta(2,:)
               zeros(1,size(zeta,2))                                                         ,... %dN3/dzeta(1,:)
              (16*zeta(2,:))./3 + (2*sqrt(3))/9                            ;    %dN3/dzeta(2,:)
              -8*zeta(1,:)                                                 ,... %dN4/dzeta(1,:)
              (8*zeta(2,:))./3 - (8*sqrt(3))/9                             ;    %dN4/dzeta(2,:)
              (8*sqrt(3)*zeta(2,:))./3 + 4/3                               ,... %dN5/dzeta(1,:)
              (8*sqrt(3)*zeta(1,:))./3 - (16*zeta(2,:))./3 + (4*sqrt(3))/9 ;    %dN5/dzeta(2,:)
             -(8*sqrt(3)*zeta(2,:))./3 - 4/3                               ,... %dN6/dzeta(1,:)
              (4*sqrt(3))/9 - (8*sqrt(3)*zeta(1,:))./3 - (16*zeta(2,:))./3 ];   %dN6/dzeta(2,:)

% loop over each element
for iel = 1:Nel
    
    % node coordinates
    xel = [node(conn(iel,:),1), node(conn(iel,:),2)];
    
    % compute distances
    R = sqrt( (repmat(xel(:,1),1,6)-repmat(xel(:,1)',6,1)).^2 + ...
              (repmat(xel(:,2),1,6)-repmat(xel(:,2)',6,1)).^2 );
    R = max(R(:));
    S = ceil(boundrange*R);
    
    % Dilation kernel
    D = ones(S,S);
    
    % find pixels inside element
    xv = node(conn(iel,[1 4 2 5 3 6]),1);
    yv = node(conn(iel,[1 4 2 5 3 6]),2);
    In = inpolygon(X,Y,xv,yv);
    In = conv2(double(In),D,'same');
    Ipx = find(In>0);
    Npx = length(Ipx);
    xpx = [X(Ipx) Y(Ipx)];
    
    % initial guess
    % ==============
    
    % create a space in master coordinates (bigger than element)
    z1 = linspace(-0.7,0.7,10*S);
    z2 = linspace(-0.5,0.8,10*S);
    [Z1, Z2] = meshgrid(z1,z2);
    
    % compute mapped real coordinates
    Xel = N(L([Z1(:)';Z2(:)']))'*xel;
    
    % interpolated master coordinates in pixel locations
    Z1 = griddata(Xel(:,1),Xel(:,2),Z1(:),xpx(:,1),xpx(:,2),'cubic');
    Z2 = griddata(Xel(:,1),Xel(:,2),Z2(:),xpx(:,1),xpx(:,2),'cubic');
    z = [Z1(:)'; Z2(:)'];

    % storing for comparison
    zinit = z;
    finit = (N(L(z))' * xel - xpx);
    
    % Iterative Solving
    % ==============
    % vectorized (all pixels at once)
    if itsolve
        for i = 1:maxit
            f = N(L(z))' * xel - xpx;
            df = dN(z)'*xel;
            df = sparse([1:2*Npx, 1:2*Npx],[1:Npx, 1:Npx, Npx+(1:Npx), Npx+(1:Npx)],df(:),2*Npx,2*Npx);
            dz = df\-f(:);
            dz = reshape(dz,Npx,2)';
            z = z + dz;
            if norm(f) < convcrit
                break
            end
        end
        
        if i == maxit
            % check where the initial guess was better
            Iinit = find( abs(f) > abs(finit) );
            z(Iinit) = zinit(Iinit);
            f = N(L(z))' * xel - xpx;
            
            fprintf('max iterations, el:%d, max err:%10.3e [px]\n',iel,max(abs(f(:))));
        end
    end
    
    % find pixels which were inside an element
    l = L(z);
    Iin = find(l(1,:)>=-tol & l(1,:)<=1+tol &...
               l(2,:)>=-tol & l(2,:)<=1+tol &...
               l(3,:)>=-tol & l(3,:)<=1+tol);
    % add T6 values to the basis matrix
    phi(Ipx(Iin),conn(iel,:)) = N(L(z(:,Iin)))';

    % update the waitbar
    if mod(iel,floor(Nel/10)) == 1
        bcwaitbar(H,iel/Nel,sprintf('building T6 basis (%d/%d)',iel,Nel));
    end

end
bcwaitbar(H);
