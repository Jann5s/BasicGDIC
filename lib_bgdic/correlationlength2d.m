function [zeta, theta, A] = correlationlength2d(f,varargin)
%CORRELATIONLENGTH2D calculates the characteristic length within which data
%    points can be considered correlated. Internally first the 2D
%    autocorrelation function is computed using the fast fourier transfer
%    [1], which is then interpolated on a cylindrical coordinate system and
%    then numerically inverted for each angle theta to find the length at
%    which the ACF intersects 1/exp(1), which is the correlation length.
%
%    The image f is assumed to be defined at integer discretization steps,
%    i.e. x = 1:m, y = 1:n, where [n, m] = size(f).
%
%    references:
%    [1] http://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
%
%    Default:
%    zeta = CORRELATIONLENGTH2D(f)
%    uses 1/exp(1) as the threshold/intersection of the ACF
%    returns the average correlation length (averaged over theta)
%
%    Optional Input:
%    zeta = CORRELATIONLENGTH2D(f,th)
%    uses th as the threshold/intersection of the ACF
%
%    Optional Output:
%    [Z, T] = CORRELATIONLENGTH2D(f)
%    returns the correlation lengths (Z) at each angle (T)   
%
%    See also FFT2, IFFT2

% one or two input variables
if nargin == 1
    th = 1/exp(1);
else
    th = varargin{1};
end

% convert to grayvalue (and to double precision)
if ndims(f) == 3
    f = 0.2989 * double(f(:,:,1)) + 0.5870 * double(f(:,:,2)) + 0.1140 * double(f(:,:,3));
else
    f = double(f);
end

% periodic mirroring (to reduce antialiassing)
if true
    f = [f fliplr(f) ; flipud(f) rot90(f,2)];
end

% number of datapoints
[n, m] = size(f);

% correct for the mean of the pattern
f = f - mean(f(:));

% signal standard deviation (note, not the rms)
sigma = std(f(:));

% fast fourier transform of f
F = fft2(f);

% fast calculation of the autocorrelation function (correctly normalized)
acf = ifft2( conj(F) .* F ) ./ (n * m * sigma^2) ;

% shifted acf
acfs = fftshift(acf);

% acf coordinate space (weird because the way fft returns the results)
ax = (1:m);
ay = (1:n);
ax = ax - ax(floor(m/2)+1);
ay = ay - ay(floor(n/2)+1);
[Ax, Ay] = meshgrid(ax,ay);

% cutoff for the ACF peak (should be slighly bigger than zero)
cutoff = 0.9*th;

% find and number islands of values bigger than cutoff
B = bwlabel(acfs >= cutoff);
% find the center island 
mnx = round(m/2);
mny = round(n/2);
centind = B(mny,mnx);
B = (B==centind);

% get the center island indices
I = find(B==1);

% get the acf peak points
x = Ax(I);
y = Ay(I);

% Circular interpolation
% =============================

% maximum circle coordinate
rmax = max(sqrt(x.^2 + y.^2));

% number of coordinates
Ntheta = 60;
Nr = min([10 3*ceil(rmax)]);

% create cylindrical space
r = linspace(0,rmax,Nr);
theta = linspace(0,2*pi,2*Ntheta);
[T, R] = meshgrid(theta,r);

% and cartesian expression of cylindrical space
X = R.*cos(T);
Y = R.*sin(T);

% Interpolate the acf on the cylindrical space
Z = interp2(Ax,Ay,acfs,X,Y,'spline');

% average values at pi angle (acf is symetrical)
if true
    Z = 0.5*(Z(:,1:Ntheta)+Z(:,Ntheta+(1:Ntheta)));
    R = R(:,1:Ntheta);
    theta = theta(1:Ntheta);
else
    % do not use symetry
    Ntheta = 2*Ntheta;
end

% get the correlation lenght (invert the radial acf profiles)
zeta = zeros(Ntheta,1);
for k = 1:Ntheta
    diff(Z(:,k));
    I = find(diff(Z(:,k))>0,1);
    if isempty(I)
        I = Nr;
    end
    zeta(k) = interp1(Z(1:I,k),R(1:I,k),th,'spline','extrap');
end


% one or two output variables
if nargout == 1
    zeta = mean(zeta);
elseif nargout == 2;
    zeta = zeta(:);
    theta = theta(:);
elseif nargout == 3;
    zeta = zeta(:);
    theta = theta(:);
    A.A = acfs;
    A.x = ax;
    A.y = ay;
end


% Plotting 
% =================================
if nargin <= 2
    return
end

colors = jncurvemap(Ntheta);

figure;
imagesc(f)
daspect([1 1 1]);
colorbar
title('image (mirrored to reduce antialiassing)')
xlabel('x [px]')
ylabel('y [px]')

figure; 
plot3(Ax(I),Ay(I),acfs(I),'.')
hold on
h = plot3(X(:,1:Ntheta),Y(:,1:Ntheta),Z);
set(h,{'Color'},colors);
daspect([1 1 1]);
title('autocorrelation function (selected data points) with radial interpolation')
xlabel('ux [px]')
ylabel('uy [px]')

figure;
h = plot(R,Z,'.-');
hold on
hz = plot(zeta,th,'s');
plot((1:m)-1,acf(1,:),'o-',(1:n)-1,acf(:,1),'d-')
set(h,{'Color'},colors);
set(hz,{'Color'},colors);
set(gca,'xlim',[0,2*rmax])
title('autocorrelation function for various angles theta')
xlabel('R [px]')
ylabel('ACF')

figure;
imagesc(ax,ay,acfs)
daspect([1 1 1]);
colorbar
set(gca,'xlim',2*rmax*[-1,1])
set(gca,'ylim',2*rmax*[-1,1])
hold on
plot([zeta zeta].*cos([theta theta+pi]),[zeta zeta].*sin([theta theta+pi]),'-k')
title('autocorrelation function (with correlation contour)')
xlabel('ux [px]')
ylabel('uy [px]')

figure;
h = plot(theta*180/pi,zeta,'.-');
title('autocorrelation length for various angles theta')
xlabel('\theta [deg]')
ylabel('\zeta')


