function [A x y] = coarsegrain(A,varargin)
%COARSEGRAIN creates a coarse-grained (superpixeled) version of an image 
%    [A x y] = coarsegrain(A,x,y,level,options)
%
%    INPUTS:
%    A = (n x m) image matrix 
%        RGB is converted to grayscale
%        class is converted to double
%
%    OPTIONAL INPUTS:
%    x = 1:m
%    y = 1:n
%    level = 1;
%    maskthreshold = []; 
%        maskthreshold requires a level to also be specified
%        superpixels > maskthreshold = true
%        a zero threshold makes the mask grow when coarse graining
%        a one threshold makes the mask shrink when coarse graining
%
%    the image matrix must always be the first input
%
%    EXAMPLE:
%    f = coarsegrain(f,x,y,3);
%
%    EXAMPLE:
%    Im = maskimage(mask,n,m);
%    Im = coarsegrain(Im,3);
%
%    See also MASKIMAGE, MASKINDEX, MASKINVERT, MASKPLOT.
%
% copyright: Jan Neggers, 2013


% Function warning message id (for suppressing warnings)
% ======================
msgid = 'JNDIC:coarsegrain';

% input handling
% ======================

% get the size of the image
[n m c] = size(A);

% initiate the options
x = [];
y = [];
level = [];
maskthreshold = [];

% read the varargin cell array
for k = 1:length(varargin)
    if length(varargin{k}) == m
        x = varargin{k};
    elseif length(varargin{k}) == n
        y = varargin{k};
    elseif isempty(level) && length(varargin{k}) == 1
        level = varargin{k};
    elseif ~isempty(level) && length(varargin{k}) == 1
        maskthreshold = varargin{k};
    end
end

% fill in the defaults
if isempty(x)
    x = 1:m;
end
if isempty(y)
    y = 1:n;
end
if isempty(level)
    level = 1;
end
if isempty(maskthreshold) && islogical(A)
    maskthreshold = 0;
end

% force x and y to be rows
x = x(:).';
y = y(:).';

% convert RGB to grayscale
if c == 3
    warning(msgid,'coarsegrain: converting RGB image to grayscale')
    A = rgb2gray(A);
end

% cg-level = 0 ==> nothing todo
if level == 0
    return
end

% set image to double
A = double(A);

% any NaN's in the images
nantest = any(isnan(A(:)));

% Cropping the image to integer superpixels
% ========================

% superpixelsize
sps = 2^level;

% number of coarse grain pixels (after last iteration)
cgn = floor(n/sps);
cgm = floor(m/sps);

% respective number of pixels in fine grain image
fgn = cgn * sps;
fgm = cgm * sps;

% number of crop pixels
cropn = n - fgn;
cropm = m - fgm;

% split crop in left and right part (and top and bottom)
cropy(1) = cropn - round(cropn/2);
cropy(2) = cropn - cropy(1);
cropx(1) = cropm - round(cropm/2);
cropx(2) = cropm - cropx(1);

% crop the image
In = cropy(1)+1:n-cropy(2);
Im = cropx(1)+1:m-cropx(2);
A = A(In,Im);
[n m] = size(A);

% crop x and y
x = x(Im);
y = y(In);


if nantest
    % coarse graining with NaN's (slower)
    % ========================
    for k = 1:level
        % number of coarse grain pixels
        cgn = n/2;
        cgm = m/2;

        % initialize the image matrices
        cgA = zeros(cgn,cgm);
        cgx = zeros(1,cgm);
        cgy = zeros(1,cgn);
        
        % initialize matrices which count non-NaN's
        sumA = zeros(cgn,cgm);
        
        for in = 1:2
            for im = 1:2
                % get the same subpixel in each superpixel
                In = in:2:n;
                Im = im:2:m;
                
                % get all the subpixels in one matrix
                Asub = A(In,Im);
                % only add x over the rows
                if in == 1;  xsub = x(Im); else xsub = 0; end
                % only add y over the columns
                if im == 1;  ysub = y(In); else ysub = 0; end
                
                % add that subpixel to the previous pixel (only if it is not NaN)
                cgA(~isnan(Asub)) = cgA(~isnan(Asub)) + Asub(~isnan(Asub));
                cgx = cgx + xsub;
                cgy = cgy + ysub;
                
                % count the number of pixels used (non-NaN's)
                sumA  = sumA + ~isnan(Asub);
            end
        end
        
        % devide by the number of pixels (excluding NaN's)
        A = cgA ./ sumA;
        x = cgx ./ 2;
        y = cgy ./ 2;
        n = cgn;
        m = cgm;
    end
else
    % coarse graining without NaN's (faster)
    % ========================
    for k = 1:level
        % number of coarse grain pixels
        cgn = n/2;
        cgm = m/2;
        
        % construct the permutation matrices
        Pn = repmat(eye(cgn),2,1);
        Pn = reshape(Pn,cgn,n)';
        Pm = repmat(eye(cgm),2,1);
        Pm = reshape(Pm,cgm,m)';
        
        % coarse grain x and y
        x = ( x * Pm ) ./ 2;
        y = ( y * Pn ) ./ 2;
        
        % coarse grain the image
        A = ( Pn' * A * Pm ) ./ 4;
        n = cgn;
        m = cgm;
    end
end

% apply maskthreshold
% ========================
if ~isempty(maskthreshold)
    A = A > maskthreshold;
end