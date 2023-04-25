function colors = jncurvemap(N,varargin)
% colors = jncurvemap(N)
%
% creates a cell-array of colors of length N. Note, this is not the same as
% the typical matrices obtained from a colormap type command (e.g. jet).
% the colors cell-array can be directly applied to a list of handles.
%
% example:
% x = linspace(0,1,100);
% N = 20;
% a = 1:N;
% for k = 1:N
%     hold on
%     h(k) = plot(x,x.^a(k));
% end
% colors = jncurvemap(N);
% set(h,{'color'},colors)
%
% optional:
% colors = jncurvemap(N,[H1 H2],[S1 S2],[L1 L2])
% where ?1 is the start point (i.e. first color) and ?2 the end point,
% H is the hue (i.e. color):
%          H=0   -> red, 
%          H=0.3 -> green,
%          H=0.6 -> blue,
%          H=1   -> red,   
% S is the saturation (i.e. color intensity)
%          S=0   -> no color -> grayscale, 
%          S=1   -> max color -> full color intensity,
% L is the lightness (i.e. brightness)
%          S=0   -> black, 
%          S=1   -> white,
% default values are
% colors = jncurvemap(N,[0.9 0.3],[0.7 0.7],[0.2 0.4])
%
% Note, in general 0 <= H,S,L <= 1, however for H > 1 the colors are looped
% around, making it possible to go around the color circle boundaries, or
% even loop around more than once.

% defaults
H1 = 0.9;
H2 = 0.3;
S1 = 0.7;
S2 = 0.7;
L1 = 0.2;
L2 = 0.4;

if nargin == 4
    H1 = varargin{1}(1);
    H2 = varargin{1}(2);
    S1 = varargin{2}(1);
    S2 = varargin{2}(2);
    L1 = varargin{3}(1);
    L2 = varargin{3}(2);
elseif nargin == 0
    % give a small demo
    help jncurvemap
    figure;
    N = 20;
    x = linspace(0,1,50);
    a = linspace(0.5,1,N);
    for k = 1:N
        hold on
        h(k) = plot(x,a(k)*sin(pi*x),'o-','linewidth',2);
    end
    colors = jncurvemap(N);
    set(h,{'color'},colors);
    error('please specify how many colors to produce: colors = jncurvemap(N)')
end

% create color spaces
H = linspace(H1,H2,N);
S = linspace(S1,S2,N);
L = linspace(L1,L2,N);

% loop around for H > 1;
H = H + floor(min(H));
H = H - floor(H);

% create color matrix
HSL = [H(:) S(:) L(:)];

% convert to RGB
RGB = hsl2rgb(HSL);

% format the output
colors = RGB;
colors = mat2cell(colors,ones(1,N),3);

function rgb=hsl2rgb(hsl)

%Converts Hue-Saturation-Luminance Color value to Red-Green-Blue Color 
%value
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

rgb=round(rgb.*100000)./100000; 
%Sometimes the result is 1+eps instead of 1 or 0-eps instead of 0 ... so to
%get rid of this I am rounding to 5 decimal places)