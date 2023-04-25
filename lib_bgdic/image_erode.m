function A = image_erode(A,R)
% A = erode(A,R), will erode the mask A with radius R. A is expected to
% be a logical matrix (i.e. with zeros and ones)
%
% See also, erode

% a small number
eps = 1e-6;

% make sure A a double
A = double(A);

% Create a Disk shaped kernel
s = 2*ceil(R+0.5)-1;
x = (1:s)-mean(1:s);
[X, Y] = meshgrid(x,x);
K = zeros(s,s);
K( X.^2+Y.^2 <= R^2) = 1;
K = K./sum(K(:));

% Convolve with K
A = conv2(A,K,'same');

% The dilated mask is everywhere where A is untouched by the kernel 
% (i.e. A == 1)
A = abs(A - 1) < eps;