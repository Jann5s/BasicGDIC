function cmap = jnmap(varargin)
% Create Colormaps with monotomic change in Luminace and change in
% GrayValue (if converted).
%
% cmap = jnmap(N) : produces and Nx3 matrix as do the default matlab
%                   colormaps
%
% cmap = jnmap(N,colorname) : same but selects a color scheme with name
%                             colorname
%
% currently available colorschemes are
% - Blueish (b)
% - Greenish (g)
% - Redish (r)
%
% examples:
% cmap = jnmap(64,'Greenish')
% cmap = jnmap(64,'g')


if nargin == 0
    N = 64;
    cname = 'blueish';
elseif nargin == 1
    if isnumeric(varargin{1})
        N = varargin{1};
        cname = 'blueish';
    else
        N = 64;
        cname = varargin{1};
    end
else
    N = varargin{1};
    cname = varargin{2};
end

% get the colormap
if any(strcmpi(cname,{'blueish','blue','b'}))
    % Blueish
    A = 24/36;
    B = 15/36;
    C = .5;
elseif any(strcmpi(cname,{'greenish','green','g'}))
    % Greenish
    A = 19/36;
    B = 6/36;
    C = .5;
elseif any(strcmpi(cname,{'redish','red','r'}))
    % Redish
    A = 30/36;
    B = 45/36;
    C =  .5;
else
    error('unknown colormap');
end

H = linspace(A,B,N);
H = H - floor(H);
S = linspace(C,C,N);
L = linspace(0,1,N);

HSL = [H(:) S(:) L(:)];

RGB = hsl2rgb(HSL);
cmap = RGB;



function rgb=hsl2rgb(hsl)

%Converts Hue-Saturation-Luminance Color value to Red-Green-Blue Color value
%
%Usage
%       RGB = hsl2rgb(HSL)
%
%   converts HSL, a M X 3 color matrix with values between 0 and 1
%   into RGB, a M X 3 color matrix with values between 0 and 1
%
%See also rgb2hsl, rgb2hsv, hsv2rgb

%Suresh E Joel, April 26,2003

if nargin<1,
    error('Too few arguements for hsl2rgb');
    return;
elseif nargin>1,
    error('Too many arguements for hsl2rgb');
    return;
end;

if max(max(hsl))>1 | min(min(hsl))<0,
    error('HSL values have to be between 0 and 1');
    return;
end;

for i=1:size(hsl,1),
    if hsl(i,2)==0,%when sat is 0
        rgb(i,1:3)=hsl(i,3);% all values are same as luminance
    end;
    if hsl(i,3)<0.5,
        temp2=hsl(i,3)*(1+hsl(i,2));
    else
        temp2=hsl(i,3)+hsl(i,2)-hsl(i,3)*hsl(i,2);
    end;
    temp1=2*hsl(i,3)-temp2;
    temp3(1)=hsl(i,1)+1/3;
    temp3(2)=hsl(i,1);
    temp3(3)=hsl(i,1)-1/3;
    for j=1:3,
        if temp3(j)>1, 
            temp3(j)=temp3(j)-1; 
        elseif temp3(j)<0, 
            temp3(j)=temp3(j)+1; 
        end;
        if 6*temp3(j)<1,
            rgb(i,j)=temp1+(temp2-temp1)*6*temp3(j);
        elseif 2*temp3(j)<1,
            rgb(i,j)=temp2;
        elseif 3*temp3(j)<2,
            rgb(i,j)=temp1+(temp2-temp1)*(2/3-temp3(j))*6;
        else
            rgb(i,j)=temp1;
        end;
    end;
end;

rgb=round(rgb.*100000)./100000; %Sometimes the result is 1+eps instead of 1 or 0-eps instead of 0 ... so to get rid of this I am rounding to 5 decimal places)