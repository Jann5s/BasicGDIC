function phi = buildphi_zernike(x,y,maxp,roi,H,varargin)

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

Im = find(x < roi(1) | x > roi(2));
In = find(y < roi(3) | y > roi(4));

% width and height of the roi
w = roi(2)-roi(1);
h = roi(4)-roi(3);
% length of the diagonal (such that rho = 1 in the corner of the roi)
L = sqrt(w^2 + h^2);

% create a space
x = 2*(x - 0.5*(roi(1)+roi(2))) / L;
y = 2*(y - 0.5*(roi(3)+roi(4))) / L;
[X,Y] = meshgrid(x,y);
r = sqrt(X(:).^2 + Y(:).^2);
t = atan2(Y(:),X(:));
Nr = length(r);

% max order
maxp = max([maxp 0]);
% number of polynomials
Nz = sum(0:maxp+1);
% normalize (such that they integrate to 1)
normalize = false;

% list of polynomials with n and m
p = (0:Nz-1)';
n = ceil((-3+sqrt(9+8*p))/2);
m = 2*p - n.*(n+2);
m_abs = abs(m);


% R polynomials
% ================
Rpow = zeros(Nr,maxp);
Rpow(:,1) = 1;
Rpow(:,2) = r;
for k = 2:maxp
    Rpow(:,k+1) = r.^k;
end

% Z polynomials
% ================
phi = zeros(Nr,Nz);
Znorm = sqrt((1+(m_abs>0)).*(n+1)./pi);
for kn = 1:Nz
    S = 0:(n(kn)-m_abs(kn))/2;
    pows = n(kn):-2:m_abs(kn);
    for kp = length(S):-1:1
        P = (1-2*mod(S(kp),2))*                ...
            prod(2:(n(kn)-S(kp)))/              ...
            prod(2:S(kp))/                     ...
            prod(2:((n(kn)-m_abs(kn))/2-S(kp)))/ ...
            prod(2:((n(kn)+m_abs(kn))/2-S(kp)));
        phi(:,kn) = phi(:,kn) + P*Rpow(:,pows(kp)+1);
    end
    
    if normalize
        phi(:,kn) = phi(:,kn)*Znorm(kn);
    end
    
    bcwaitbar(H,kn/Nz,sprintf('building zernike basis (%d/%d)',kn,Nz));
end

% Zernike polynomials
% ================
phi(:,m>0) = phi(:,m>0).*cos(t*m_abs(m>0)');
phi(:,m<0) = phi(:,m<0).*sin(t*m_abs(m<0)');

if strcmpi(precision,'single')
    phi = single(phi);
end

bcwaitbar(H);
