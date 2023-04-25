function phi = buildphi_pseudozernike(x,y,maxp,roi,H,varargin)

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
% normalize (such that they integrate to 1)
normalize = false;
% scale to -1,1
scale = false;

if scale
%     Ic = find(r <= 1);
    roi([1,2]) = 2*(roi([1,2]) - 0.5*(roi(1)+roi(2))) / L;
    roi([3,4]) = 2*(roi([3,4]) - 0.5*(roi(3)+roi(4))) / L;
    Ic = find(X > roi(1) & X < roi(2) & Y > roi(3) & Y < roi(4));
end

% number of polynomials
Nz = sum(1:2:2*(maxp+1));
p = (1:Nz)' - 1;

% indices
n = floor(sqrt(p));
m = p - n.^2 - n;
m_abs = abs(m);

% pseudo Z polynomials
% ================
phi = zeros(Nr,Nz);
for k = 1:Nz
    for s = 0:(n(k)-m_abs(k))
        D = (-1)^(s) * factorial(2*n(k) + 1 - s) ./ ...
            (factorial(s) * factorial(n(k) - m_abs(k) - s) .* ...
            factorial(n(k) + m_abs(k) - s + 1));
        
        phi(:,k) = phi(:,k) + D*r.^(n(k) - s);
    end
    
    if m(k) > 0
        phi(:,k) = phi(:,k) .* cos(m_abs(k)*t);
    elseif m(k) < 0
        phi(:,k) = phi(:,k) .* sin(m_abs(k)*t);
    end
    
    if normalize
        phi(:,k) = phi(:,k) * 2*(n(k) + 1)/pi;
    end
    
    if scale && (k > 1)
        phi(:,k) = (phi(:,k) - min(phi(Ic,k))) ./ (max(phi(Ic,k)) - min(phi(Ic,k)));
        phi(:,k) = 2*phi(:,k) - 1;
    end
    
    bcwaitbar(H,k/Nz,sprintf('building pseudo zernike basis (%d/%d)',k,Nz));
end

if strcmpi(precision,'single')
    phi = single(phi);
end

bcwaitbar(H);
