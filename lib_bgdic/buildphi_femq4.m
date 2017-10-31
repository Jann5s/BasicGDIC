function phi = buildphi_femq4(x,y,node,conn,H,varargin)
% phi = buildphi_femq4(x,y,node,conn)
%
% Function for creating computing the basis function values for each pixel
% for each degree of freedom (node), for a mesh with 4 noded quadrilateral
% elements (aka Q4).
%
% x (1,m) and y (1,n) are the pixel positions
% node is a matrix with nodal coordinates (Ndof,2)
% conn is a matrix with element connectivity (Nel,4)
%
% the node numbering convention used is
%     4 ---------- 3
%     |            |
%     |            |
%     |            |
%     1 ---------- 2

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

convcrit = 1e-10;
maxit = 40;
init = [0 ; 0];
tol = 1e-6; % pixel position tolerance

% Q4 shape functions
N = @(zeta) [0.25*(1-zeta(1,:)).*(1-zeta(2,:));
             0.25*(1+zeta(1,:)).*(1-zeta(2,:));
             0.25*(1+zeta(1,:)).*(1+zeta(2,:));
             0.25*(1-zeta(1,:)).*(1+zeta(2,:))];

% Q4 derivatives
dN = @(zeta) [ 0.25*( -1 + zeta(2,:) ), 0.25*( -1 + zeta(1,:) );
               0.25*(  1 - zeta(2,:) ), 0.25*( -1 - zeta(1,:) );
               0.25*(  1 + zeta(2,:) ), 0.25*(  1 + zeta(1,:) );
               0.25*( -1 - zeta(2,:) ), 0.25*(  1 - zeta(1,:) )];

% loop over each element
for iel = 1:Nel
    
    % node coordinates
    xel = [node(conn(iel,:),1), node(conn(iel,:),2)];

    % compute element size
    R = sqrt( (repmat(xel(:,1),1,4)-repmat(xel(:,1)',4,1)).^2 + ...
              (repmat(xel(:,2),1,4)-repmat(xel(:,2)',4,1)).^2 );
    R = ceil(max(R(:)));
    
    % find pixels inside element
    xv = node(conn(iel,:),1);
    yv = node(conn(iel,:),2);
    Ipx = inpolygon(X,Y,xv,yv);
    Ipx = find(Ipx);
    Npx = length(Ipx);
    xpx = [X(Ipx) Y(Ipx)];

    % initial guess
    % ==============
    
    if false
        % zero initial guess (boring)
        z = zeros(2,Npx);
        
        % storing for comparison
        zinit = z;
        Finit = (N(z)' * xel - xpx);
        
    else
        % interpolation initial guess (fun)
        % create a space in master coordinates (bigger than element)
        z1 = linspace(-1.1,1.1,round(2*R));
        z2 = linspace(-1.1,1.1,round(2*R));
        [Z1, Z2] = meshgrid(z1,z2);
        
        % compute mapped real coordinates
        Xel = N([Z1(:)';Z2(:)'])'*xel;
        
        % interpolated master coordinates in pixel locations
        Z1 = griddata(Xel(:,1),Xel(:,2),Z1(:),xpx(:,1),xpx(:,2),'linear');
        Z2 = griddata(Xel(:,1),Xel(:,2),Z2(:),xpx(:,1),xpx(:,2),'linear');
        z = [Z1(:)'; Z2(:)'];
        
        % storing for comparison
        zinit = z;
        finit = (N(z)' * xel - xpx);
        
    end

    % iterative solving
    % ==============
    for i = 1:maxit
        f = N(z)' * xel - xpx;
        df = dN(z)'*xel;
        df = sparse([1:2*Npx,1:2*Npx],[1:Npx,1:Npx,Npx+(1:Npx),Npx+(1:Npx)],df(:),2*Npx,2*Npx);
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
        f = N(z)' * xel - xpx;

        fprintf('max iterations, el:%d, max err:%10.3e [px]\n',iel,max(abs(f(:))));
    end
    
    % find pixels which were inside an element
    Iin = find(abs(z(1,:))<=1+tol & abs(z(2,:))<=1+tol);
    % add Q4 values to the basis matrix
    phi(Ipx(Iin),conn(iel,:)) = N(z(:,Iin))';

    % update the waitbar
    if mod(iel,floor(Nel/10)) == 1
        bcwaitbar(H,iel/Nel,sprintf('building Q4 basis (%d/%d)',iel,Nel));
    end

end
bcwaitbar(H);
