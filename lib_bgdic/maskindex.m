function I = maskindex(Im)
%MASKINDEX converts a logical image of size (n x m) to a list of indices
%
%    I = maskindex(Im)
%
%    See also MASKIMAGE, MASKINVERT.

I = find(Im);