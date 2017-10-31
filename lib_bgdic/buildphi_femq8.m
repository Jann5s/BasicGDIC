function phi = buildphi_femq8(x,y,node,conn,H,varargin)
% phi = buildphi_femq8(x,y,node,conn)
%
% Function for creating computing the basis function values for each pixel
% for each degree of freedom (node), for a mesh with 8 noded quadrilateral
% elements (aka Q8).
%
% x (1,m) and y (1,n) are the pixel positions
% node is a matrix with nodal coordinates (Ndof,2)
% conn is a matrix with element connectivity (Nel,8)
%
% the node numbering convention used is
%     4 --- 7 --- 3
%     |           |
%     8           6
%     |           |
%     1 --- 5 --- 2

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


% Q8 shape functions
N = @(z) [-0.25*(1-z(1,:)) .* (1-z(2,:)) .* (1+z(1,:)+z(2,:));
    -0.25*(1+z(1,:)) .* (1-z(2,:)) .* (1-z(1,:)+z(2,:));
    -0.25*(1+z(1,:)) .* (1+z(2,:)) .* (1-z(1,:)-z(2,:));
    -0.25*(1-z(1,:)) .* (1+z(2,:)) .* (1+z(1,:)-z(2,:));
    0.5*(1-z(1,:).^2) .* (1-z(2,:)   );
    0.5*(1+z(1,:)   ) .* (1-z(2,:).^2);
    0.5*(1-z(1,:).^2) .* (1+z(2,:)   );
    0.5*(1-z(1,:)   ) .* (1-z(2,:).^2)];

% Q8 Derivatives
dN = @(z) [0.25*( 2*z(1,:) +  z(2,:) -2*z(1,:).*z(2,:) -z(2,:).^2 ),...% dN1/dz1
    0.25*(   z(1,:) +2*z(2,:) -2*z(1,:).*z(2,:) -z(1,:).^2 );   % dN1/dz2
    0.25*( 2*z(1,:) -  z(2,:) -2*z(1,:).*z(2,:) +z(2,:).^2 ),...% dN2/dz1
    0.25*(-  z(1,:) +2*z(2,:) +2*z(1,:).*z(2,:) -z(1,:).^2 );   % dN2/dz2
    0.25*( 2*z(1,:) +  z(2,:) +2*z(1,:).*z(2,:) +z(2,:).^2 ),...% dN3/dz1
    0.25*(   z(1,:) +2*z(2,:) +2*z(1,:).*z(2,:) +z(1,:).^2 );   % dN3/dz2
    0.25*( 2*z(1,:) -  z(2,:) +2*z(1,:).*z(2,:) -z(2,:).^2 ),...% dN4/dz1
    0.25*(-  z(1,:) +2*z(2,:) -2*z(1,:).*z(2,:) +z(1,:).^2 );   % dN4/dz2
    -z(1,:) + z(1,:).*z(2,:), 0.5 * (z(1,:).^2 - 1)  ; % dN5/dz1, dN5/dz2
    -0.5 * (z(2,:).^2 - 1)  ,-z(2,:) - z(1,:).*z(2,:); % dN6/dz1, dN6/dz2
    -z(1,:) - z(1,:).*z(2,:),-0.5 * (z(1,:).^2 - 1)  ; % dN7/dz1, dN7/dz2
    0.5 * (z(2,:).^2 - 1)  ,-z(2,:) + z(1,:).*z(2,:)];% dN8/dz1, dN8/dz2

% loop over each element
for iel = 1:Nel
    
    % node coordinates
    xel = [node(conn(iel,:),1), node(conn(iel,:),2)];
    
    % compute distances
    R = sqrt( (repmat(xel(:,1),1,8)-repmat(xel(:,1)',8,1)).^2 + ...
              (repmat(xel(:,2),1,8)-repmat(xel(:,2)',8,1)).^2 );
    R = max(R(:));
    S = ceil(boundrange*R);
    
    % Dilation kernel
    D = ones(S,S);
    
    % find pixels inside element
    xv = node(conn(iel,[1 5 2 6 3 7 4 8]),1);
    yv = node(conn(iel,[1 5 2 6 3 7 4 8]),2);
    In = inpolygon(X,Y,xv,yv);
    In = conv2(double(In),D,'same');
    Ipx = find(In>0);
    Npx = length(Ipx);
    xpx = [X(Ipx) Y(Ipx)];
    
    % initial guess
    % ==============
    
    % create a space in master coordinates (bigger than element)
    z1 = linspace(-2,2,15*S);
    z2 = linspace(-2,2,15*S);
    [Z1, Z2] = meshgrid(z1,z2);
    
    % compute mapped real coordinates
    Xel = N([Z1(:)';Z2(:)'])'*xel;
    
    % interpolated master coordinates in pixel locations
    Z1 = griddata(Xel(:,1),Xel(:,2),Z1(:),xpx(:,1),xpx(:,2),'cubic');
    Z2 = griddata(Xel(:,1),Xel(:,2),Z2(:),xpx(:,1),xpx(:,2),'cubic');
    z = [Z1(:)'; Z2(:)'];

    % storing for comparison
    zinit = z;
    finit = (N(z)' * xel - xpx);
    
    % Iterative Solving
    % ==============
    % vectorized (all pixels at once)
    if itsolve
        for i = 1:maxit
            f = N(z)' * xel - xpx;
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
            f = N(z)' * xel - xpx;
            
            fprintf('max iterations, el:%d, max err:%10.3e [px]\n',iel,max(abs(f(:))));
        end
    end
    
    % find pixels which were inside an element
    Iin = find(abs(z(1,:))<=1+tol & abs(z(2,:))<=1+tol);
    % add Q8 values to the basis matrix
    phi(Ipx(Iin),conn(iel,:)) = N(z(:,Iin))';

    % update the waitbar
    if mod(iel,floor(Nel/10)) == 1
        bcwaitbar(H,iel/Nel,sprintf('building Q8 basis (%d/%d)',iel,Nel));
    end

end
bcwaitbar(H);

% ===============================
% replace this code with the vectorized code above to solve the problem for
% each pixel individually

% solve per pixel
% F = zeros(Npx,2);
% for ipx = 1:Npx;
%     for i = 1:maxit
%         f = N(z(:,ipx))' * xel - xpx(ipx,:);
%         F(ipx,:) = f;
%         df = dN(z(:,ipx))'*xel;
%         dz = df\-f(:);
%         z(:,ipx) = z(:,ipx) + dz;
%         if norm(f) < convcrit
%             break
%         end
%         if i == maxit
%             fprintf('max iterations, el:%d, px:%d, res:%10.3e\n',iel,ipx,norm(f));
%         end
%     end
% end