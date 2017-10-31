function [V, D] = eig2d(A)
% [V, D] = eig2d(A)
% compute the eigenvalues and eigen vectors of a matrix (field) of 2D tensors
%
% A = {2,2}(n,m), i.e. a 2x2 cell array of matrices, where
%    A{1,1}(:,:) is A11, or Axx
%    A{1,2}(:,:) is A12, or Axy
%    A{2,1}(:,:) is A21, or Ayx
%    A{2,2}(:,:) is A22, or Ayy
%
% the returned matrix V has the same structure, where
%    V{1,1}(:,:) is the first component of the first eigenvector
%    V{1,2}(:,:) is the second component of the first eigenvector
%    V{2,1}(:,:) is the first component of the second eigenvector
%    V{2,2}(:,:) is the second component of the second eigenvector
% and the matrix D has the structure
%    D{1}(:,:) is the first eigenvalue
%    D{2}(:,:) is the second eigenvalue
%
% [V, D] = eig2d(A,'gpu')

[n, m, k] = size(A{1,1});
cls = class(A{1,1});

V{1,1} = zeros(n,m,cls);
V{1,2} = zeros(n,m,cls);
V{2,1} = zeros(n,m,cls);
V{2,2} = zeros(n,m,cls);

D{1,1} = zeros(n,m,cls);
D{2,1} = zeros(n,m,cls);

% solving the characteristic polynomial
b = -A{1,1} -  A{2,2};
c =  A{1,1} .* A{2,2} - A{1,2}.*A{2,1};
d = b.^2 - 4*c;

% fix if near zero values are negative
d(d<0 & d > -eps) = 0;

d = sqrt(d);
D{1,1} = (-b + d)./2;
D{2,1} = (-b - d)./2;

% case one, if A{1,2} == 0 and A{2,1} == 0
I1 = A{1,2} == 0 & A{2,1} == 0;
V{1,1}(I1) = 1;
V{1,2}(I1) = 0;
V{2,1}(I1) = 0;
V{2,2}(I1) = 1;

% case two, A{1,2} ~= 0
I2 = A{1,2} ~= 0;
I2 = I2 & ~I1;
V{1,1}(I2) = A{1,2}(I2);
V{1,2}(I2) = D{1,1}(I2)-A{1,1}(I2);
V{2,1}(I2) = A{1,2}(I2);
V{2,2}(I2) = D{2,1}(I2)-A{1,1}(I2);

% case three, A{2,1} ~= 0
I3 = A{2,1} ~= 0;
I3 = I3 & ~I2 & ~I1;
V{1,1}(I3) = D{1,1}(I3)-A{2,2}(I3);
V{1,2}(I3) = A{2,1}(I3);
V{2,1}(I3) = D{2,1}(I3)-A{2,2}(I3);
V{2,2}(I3) = A{2,1}(I3);

% normalize
N1 = sqrt(V{1,1}.^2 + V{1,2}.^2);
N2 = sqrt(V{2,1}.^2 + V{2,2}.^2);
V{1,1} = V{1,1}./N1;
V{1,2} = V{1,2}./N1;
V{2,1} = V{2,1}./N2;
V{2,2} = V{2,2}./N2;

% sorting, find index where D1 > D2 and swap values
I = D{1,1} > D{2,1};
[D{1,1}(I), D{2,1}(I)] = deal(D{2,1}(I), D{1,1}(I));
[V{1,1}(I), V{2,1}(I)] = deal(V{2,1}(I), V{1,1}(I));
[V{1,2}(I), V{2,2}(I)] = deal(V{2,2}(I), V{1,2}(I));

% I1 = D{1,1} <= D{2,1};
% I2 = ~I1;
% [D{1,1}(I1), D{1,1}(I2), D{2,1}(I1), D{2,1}(I2)] = deal(D{1,1}(I1),D{2,1}(I2),D{2,1}(I1),D{1,1}(I2));
% [V{1,1}(I1), V{1,1}(I2), V{2,1}(I1), V{2,1}(I2)] = deal(V{1,1}(I1),V{2,1}(I2),V{2,1}(I1),V{1,1}(I2));
% [V{1,2}(I1), V{1,2}(I2), V{2,2}(I1), V{2,2}(I2)] = deal(V{1,2}(I1),V{2,2}(I2),V{2,2}(I1),V{1,2}(I2));

