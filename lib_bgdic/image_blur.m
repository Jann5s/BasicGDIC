function A = image_blur(A,R)
% B = image_blur(A,R) blurs the image A with a gaussian kernel of radius R
%
% Notes:
% - images are mirrored to deal with the edges.
% - internally, doubles are used, but B is returned as the same type as A.
%   When A is an integer, this may lead to intensity clipping.
%
% see also, image_grayscale

% Changelog: 
% 2015-03, JN

% convert to doubles to preserve accuracy
cls = class(A);
A = double(A);

% create the kernel
% ------------------------
kernel = 'gaussian';
if strcmpi(kernel,'gaussian')
    % box size
    s = 2*ceil(4*R+0.5)-1;
    x = (1:s)-mean(1:s);
    
    K = (1/(sqrt(2*pi)*R))*exp( -(x.^2)./(2*R^2) );
end

for i = 1:ndims(A)
    
    % mirroring the image
    % ------------------------
    [an,am,ap] = size(A);
    km         = length(K);
    
    % mirror distance
    Nm = floor(km/2);
    
    % skip a direction if it has too few pixels
    if Nm >= am
        if ndims(A) == 2
            % transpose the image and blur again
            A = permute(A,[2 1 3]);
        elseif ndims(A) == 3
            % transpose the image and blur again
            A = permute(A,[2 3 1]);
        end
        continue
    end
    
    Amir = zeros(an,am+2*Nm,ap);
    
    % image rows and columns in mirrored image
    Im = Nm+1:Nm+am;
        
    % image rows and columns in the original image
    Jm = 1:am;
    
    % mirror parts in mirrored image
    Iml = Nm:-1:1; % left columns
    Imr = 2*Nm+am:-1:Nm+1+am; % right columns
    
    % mirror parts in original image
    Jml = 2:Nm+1;
    Jmr = (am:am+Nm-1)-Nm;

    % Create an image with mirrored edges
    Amir(:,[Iml,Im,Imr],:) = A(:,[Jml,Jm,Jmr],:);
    
    % Apply the blur using convolution
    % ------------------------

    if ndims(A) == 2
        A = conv2(Amir,K,'valid');

        % transpose the image and blur again
        A = permute(A,[2 1 3]);
    elseif ndims(A) == 3
        A = convn(Amir,K,'valid');

        % transpose the image and blur again
        A = permute(A,[2 3 1]);
    end
    
    
end

% Recasting
% ------------------------

% cast back to the original type
A = cast(A,cls);


