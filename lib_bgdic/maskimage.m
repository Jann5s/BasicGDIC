function Im = maskimage(I,n,m)
%MASKIMAGE converts a list of indices to a logical image of size (n x m)
%
%    Im = maskimage(I,n,m)
%
%    See also MASKINDEX, MASKINVERT.

Im = false(n,m);
Im(I) = true;