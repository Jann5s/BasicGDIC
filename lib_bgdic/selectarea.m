function varargout = selectarea(varargin)
% A = selectarea, allows the user to interactively draw a shape in the
% current axes. The default shape is a rectangle. The shape control points
% are returned in matrix A, which for the case of the rectangle are the
% positions of two diagonally opposed corners, see more details below.
% Control points are moved by dragging them to a new position. The shape
% curve is finalized/confirmed by double clicking anywhere in the figure.
% Right clicking anywhere in the figure will present a context menu where
% for instance the line color can be changed.
%
% [A, P] = selectarea, also return a position matrix: P = [x(:), y(:)]
%
% A = selectarea(B), initializes the rectangle using the control points
% defined in B
%
% A = selectarea(B,shape), where shape is any of:
%   'Point'     : draw a point. B = [x,y];
%   'Square'    : draw a square from the center point and a second point 
%                 on the square. B = [xc,yc;x2,y2]; 
%   'Rectangle' : draw a rectangle using two corners as control points,
%                 B = [x1,y1;x2,y2]; 
%   'Circle'    : draw a circle from the center point and a second point 
%                 on the circle. B = [xc,yc;x2,y2];
%   'Ellipse'   : draw an ellipse from the center point and a two points 
%                 on the ellipse. The first point defines the major radius
%                 and the angle of the central ellipse axes, the second
%                 point only defines the minor radius. 
%                 B = [xc,yc;x2,y2;x3,y3];
%   'Polyline'  : draw a polyline, if the first and the last control points 
%                 coinside the polyline is closed (i.e. its a polygon)
%                 B = [x1,y1;...;xn,yn];
%   'BSpline'   : draw a bspline of order 3, if the first and the last 
%                 control points coinside the bspline is closed. 
%                 B = [x1,y1;...;xn,yn];
%
% A = selectarea(B,'BSpline',p), draw a bspline of order p
%
% A = selectarea(B,'BSpline',phi), draw a bspline using the shape function
% matrix phi (i.e. P = C*phi). phi should be of size = [N, NB], Here N is
% the number of points on the bspline and NB the number of control points.

% the the current handles
ha = gca;
hf = gcf;

% raise the current figure
figure(hf);
set(hf,'CurrentAxes',ha);

% define a list of line colors
linecol = {'r','Red';'g','Green';'b','Blue';'c','Cyan';'m','Magenta';'y','Yellow';'k','Black';'w','White'};
linewidth = [2,0.5];
markersize = 12;

% number of points for most shapes
UserData.N = 64;

% default bspline order
p = 3;

% initiate the arguments
xlim = get(ha,'xlim');
ylim = get(ha,'ylim');
rng(1) = diff(xlim);
rng(2) = diff(ylim);
if nargin == 0
    B = [xlim(1)+0.1*rng(1), ylim(1)+0.1*rng(2); xlim(1)+0.9*rng(1), ylim(1)+0.9*rng(2)];
    shape = 'rectangle';
elseif nargin == 1
    B = varargin{1};
    shape = 'rectangle';
elseif nargin == 2
    B = varargin{1};
    shape = varargin{2};
elseif nargin == 3
    B = varargin{1};
    shape = varargin{2};
    p = varargin{3};
end

% get the control point position data
UserData.Closed = false;
if strcmpi(shape,'square')
    shapefun = @ShapeSquare;
    B = B(1:2,:);
elseif strcmpi(shape,'rectangle')
    shapefun = @ShapeRectangle;
    B = B(1:2,:);
    B = [mean(B) ; B];    
elseif strcmpi(shape,'circle')
    shapefun = @ShapeCircle;
    B = B(1:2,:);
elseif strcmpi(shape,'ellipse')
    shapefun = @ShapeEllipse;
    B = B(1:3,:);
elseif strcmpi(shape,'polyline')
    shapefun = @ShapePolyLine;
    if hypot(B(end,1)-B(1,1),B(end,2)-B(1,2)) < mean(rng)*1e-6
        UserData.Closed = true;
        B = B(1:end-1,:);
    else
        UserData.Closed = false;
    end
    % add the central rigid body motion point
    B = [mean(B) ; B];    
elseif strcmpi(shape,'bspline')
    shapefun = @ShapeBSpline;
    if numel(p) == 1
        % create the bspline shape functions
        Nb = size(B,1);
        if Nb < (p+1)
            error('Too few control points (%d) for p=%d order b-spline',Nb,p)
        end
        
        t = linspace(0,1,10*Nb-1);
        if hypot(B(end,1)-B(1,1),B(end,2)-B(1,2)) < mean(rng)*1e-6
            % Closed curve
            UserData.Closed = true;
            knots = linspace(0,1,Nb);
            dk = diff(knots(1:2));
            knots = [-(fliplr(dk:dk:p*dk)), knots, 1+(dk:dk:p*dk)];
            phi = bsplines(t,knots,p);
            for k = 1:p
                phi(:,k) = phi(:,k)+phi(:,end-(p-k));
            end
            phi = phi(:,1:Nb-1);
            B = B(1:end-1,:);
        else
            % Open Curve
            UserData.Closed = false;
            Nk = Nb - p + 1;
            knots = [zeros(1,p), linspace(0,1,Nk), ones(1,p)] ;
            phi = bsplines(t,knots,p);
        end
        UserData.phi = phi;
        UserData.phiprovided = false;
    else
        % shape functions are provided
        UserData.phi = p;
        UserData.phiprovided = true;
    end
    if size(B,1) ~= size(UserData.phi,2)
        error('the number of control points (%d) does not match the basis (phi)')
    end
    % add the central rigid body motion point
    B = [mean(B) ; B];    
elseif strcmpi(shape,'point')
    shapefun = @ShapePoint;
else
    shapefun = @ShapeRectangle;
    B = B(1:2,:);
    % add the central rigid body motion point
    B = [mean(B) ; B];    
end

% Store the current figure properties
UserData.Pointer = get(gcf,'Pointer');
UserData.NextPlot = get(ha,'NextPlot');

% if an image is used, it normally prevents the context menu, this fixed
% that
hi = findobj(hf,'Type','image');
set(hi,'HitTest','off');

% set some figure properties
set(hf,'Pointer','crosshair');
set(ha,'NextPlot','Add')
set(ha,'Clipping','off'); 
set(hi,'Clipping','off'); 

% Plot the shape for the first time
% ---------------------------------------

% Compute the shape coordinates
[x, y, B] = shapefun(B,UserData);

set(0,'CurrentFigure',hf);
set(hf,'CurrentAxes',ha);


% the background line
h1 = plot(x,y,'-','linewidth',linewidth(1));
% the foreground line
h2 = plot(x,y,'-','linewidth',linewidth(2));
% the control points
if UserData.Closed
    h3 = plot(B([2:end, 2],1),B([2:end, 2],2),'s');
else
    h3 = plot(B(2:end,1),B(2:end,2),'s');
end
if strcmpi(shape,'bspline')
    set(h3,'LineStyle','--')
end
% the rigid body motion point
h4 = plot(B(1,1),B(1,2),'o');
set([h3,h4],'MarkerSize',markersize);
set([h1;h2;h3;h4],'Clipping','off')
changecolor([],[],h1,h2,h3,h4,linecol(6,:))

% Store data which is used by the subfunctions
UserData.B = B;
UserData.h1 = h1;
UserData.h2 = h2;
UserData.h3 = h3;
UserData.h4 = h4;
UserData.ha = ha;
UserData.hf = hf;
UserData.shape = shape;
UserData.shapefun = shapefun;
UserData.Rsel = 0.02*mean(rng);
UserData.lim = [xlim; ylim];
UserData.CurrentHandle = [];
UserData.p = p;

% setup the figure such that it will allow user control
set(gcf,'UserData',UserData,...
        'WindowButtonMotionFcn',@freemove,...
        'WindowButtonDownFcn',@selecthandle,...
        'WindowButtonUpFcn',@confirmhandle,...
        'KeyPressFcn',@keyPressFcn,...
        'DoubleBuffer','on');

% define the context menu
c = uicontextmenu;
set([ha;hf;h1;h2;h3;h4],'UIContextMenu',c);

% Create top-level menu item
uimenu('Parent',c,'Label','Confirm Selection','Callback',{@confirmselect,gcf});
m = uimenu('Parent',c,'Label','Line Color');
for k = 1:8
    uimenu('Parent',m,'Label',linecol{k,2},'Callback',{@changecolor,h1,h2,h3,h4,linecol(k,:)});
end
% add control points option
if strcmpi(shape,'polyline') || strcmpi(shape,'bspline') && ~UserData.phiprovided
    uimenu('Parent',c,'Label','Add Control Point','Callback',{@addcontrolpoint,gcf});
    uimenu('Parent',c,'Label','Del Control Point','Callback',{@delcontrolpoint,gcf});
    uimenu('Parent',c,'Label','Open/Close Curve','Callback',{@openclosecurve,gcf});
end

% Wait until the user has finished
% -----------------------------------------------------
waitfor(h4)

% Finalize after the user has finished
% -----------------------------------------------------
if ~ishandle(hf)
    % if the user closed the figure
    varargout{1} = [];
    varargout{2} = [];
    return
end

% cleanup/restore figure properties
set(ha,'NextPlot',UserData.NextPlot);
set(hf,'Pointer',UserData.Pointer)
set(hf,'WindowButtonMotionFcn','','WindowButtonDownFcn','','WindowButtonUpFcn','','DoubleBuffer','off');
set([ha,hf],'UIContextMenu','');
set(hi,'HitTest','on');
set(ha,'Clipping','on'); 
set(hi,'Clipping','on'); 


% get the latest state of the curve
UserData = get(gcf,'UserData');
B = UserData.B;
if nargout == 2
    [x, y, B] = shapefun(B,UserData);
end

% Remove the extra center point from the control points list (if required)
if strcmpi(shape,'rectangle')
    B = B(2:3,:);
elseif strcmpi(shape,'polyline')
    B = B(2:end,:);
elseif strcmpi(shape,'bspline')
    B = B(2:end,:);
end

% Get the function outputs
if nargout >= 1
    varargout{1} = B;
end
if nargout == 2
    varargout{2} = [x(:), y(:)];
end

% =====================================================
function [x, y, B] = ShapeSquare(B,UserData)
% distance to the center
R = hypot(B(2,2)-B(1,2),B(2,1)-B(1,1));

% angle (with horizon)
a = atan2(B(2,2)-B(1,2),B(2,1)-B(1,1));

% make a rotated square
Ix = [-1,  1, 1, -1, -1];
Iy = [-1, -1, 1,  1, -1];
x = B(1,1) + R*Ix.*cos(a) - R*Iy.*sin(a);
y = B(1,2) + R*Ix.*sin(a) + R*Iy.*cos(a);

% fix the control points
B(2,1) = B(1,1) + R*cos(a);
B(2,2) = B(1,2) + R*sin(a);
    
% =====================================================
function [x, y, B] = ShapeRectangle(B,UserData)
x = B([2,3,3,2,2],1);
y = B([2,2,3,3,2],2);
B(1,:) = mean(B(2:3,:));

% =====================================================
function [x, y, B] = ShapeCircle(B,UserData)
a = linspace(0,2*pi,UserData.N);
R = hypot(B(1,1)-B(2,1),B(1,2)-B(2,2));
x = B(1,1) + R.*cos(a);
y = B(1,2) + R.*sin(a);

% =====================================================
function [x, y, B] = ShapeEllipse(B,UserData)
a = hypot(B(2,1)-B(1,1),B(2,2)-B(1,2));
b = hypot(B(3,1)-B(1,1),B(3,2)-B(1,2));
phi = atan2(B(2,2)-B(1,2),B(2,1)-B(1,1)) ; % ellipse angle
t = linspace(0,2*pi,UserData.N);
x = B(1,1) + a.*cos(t).*cos(phi) - b.*sin(t).*sin(phi);
y = B(1,2) + a.*cos(t).*sin(phi) + b.*sin(t).*cos(phi);
B(2,1) = B(1,1) + a*cos(phi);
B(2,2) = B(1,2) + a*sin(phi);
B(3,1) = B(1,1) + b*cos(phi+0.5*pi);
B(3,2) = B(1,2) + b*sin(phi+0.5*pi);

% =====================================================
function [x, y, B] = ShapePolyLine(B,UserData)
if UserData.Closed
    x = B([2:end, 2],1);
    y = B([2:end ,2],2);
else
    x = B(2:end,1);
    y = B(2:end,2);
end
B(1,:) = mean(B(2:end,:));

% =====================================================
function [x, y, B] = ShapeBSpline(B,UserData)
x = UserData.phi*B(2:end,1);
y = UserData.phi*B(2:end,2);
B(1,:) = mean(B(2:end,:));

% =====================================================
function [x, y, B] = ShapePoint(B,UserData)
x = B(1,1);
y = B(1,2);

% =====================================================
function openclosecurve(hObject, varargin)
H = varargin{2};
UserData = get(H,'UserData');
B = UserData.B;
B = B(2:end,:);

if UserData.Closed == true
    UserData.Closed = false;
    if strcmpi(UserData.shape,'bspline')
        B = [B ; B(1,:)];
    end
else
    UserData.Closed = true;
end

if strcmpi(UserData.shape,'bspline')
    p = UserData.p;
    Nb = size(B,1);
    t = linspace(0,1,10*Nb-1);
    if UserData.Closed
        % Closed curve
        knots = linspace(0,1,Nb);
        dk = diff(knots(1:2));
        knots = [-(fliplr(dk:dk:p*dk)), knots, 1+(dk:dk:p*dk)];
        phi = bsplines(t,knots,p);
        for k = 1:p
            phi(:,k) = phi(:,k)+phi(:,end-(p-k));
        end
        phi = phi(:,1:Nb-1);
        B = B(1:end-1,:);
    else
        % Open Curve
        Nk = Nb - p + 1;
        knots = [zeros(1,p), linspace(0,1,Nk), ones(1,p)] ;
        phi = bsplines(t,knots,p);
    end
    UserData.phi = phi;
end
B = [mean(B) ; B];

% update the curve coordinates
[x, y, B] = UserData.shapefun(B,UserData);

% update the userdata for possible changes in B
UserData.B = B;
set(H,'UserData',UserData);

% change the coordinates of the plotted curves
set(UserData.h1,'XData',x,'YData',y);
set(UserData.h2,'XData',x,'YData',y);
if UserData.Closed
    set(UserData.h3,'XData',B([2:end, 2],1),'YData',B([2:end, 2],2));
else
    set(UserData.h3,'XData',B(2:end,1),'YData',B(2:end,2));
end
set(UserData.h4,'XData',B(1,1),'YData',B(1,2));

% =====================================================
function addcontrolpoint(hObject, varargin)
H = varargin{2};
UserData = get(H,'UserData');

% get the current location
p = get(UserData.ha,'CurrentPoint');
p = p(1,1:2);

B = UserData.B;
B = B(2:end,:);
N = size(B,1);

% distance from the mouse to each control point
D = hypot(B(1:end,1)-p(1),B(1:end,2)-p(2));
% find the two closest control points
[D, I] = sort(D);

if UserData.Closed == true
    % was the mouse closer to the next point or the previous point
    if I(1) == 1
        % the previous and next node
        Ip = [N, 1+1];
        
        % compare the distance to the line
        L1 = B([I(1) Ip(1)],:);
        d1 = linepointdistance(L1,p);
        L2 = B([I(1) Ip(2)],:);
        d2 = linepointdistance(L2,p);
        
        % pick the one with the minimum difference
        if d1 < d2
            B = [B(1:N,:) ; p];
        else
            B = [B(1,:) ; p ; B(2:end,:)];
        end
    elseif I(1) == N
        % the previous and next node
        Ip = [N-1, 1];
        
        % compare the distance to the line
        L1 = B([I(1) Ip(1)],:);
        d1 = linepointdistance(L1,p);
        L2 = B([I(1) Ip(2)],:);
        d2 = linepointdistance(L2,p);
        
        % pick the one with the minimum difference
        if d1 < d2
            B = [B(1:N-1,:) ; p ; B(N,:)];
        else
            B = [B(1:N,:) ; p];
        end
    else
        Ip = [I(1)-1, I(1)+1];
        
        % compare the distance to the line
        L1 = B([I(1) Ip(1)],:);
        d1 = linepointdistance(L1,p);
        L2 = B([I(1) Ip(2)],:);
        d2 = linepointdistance(L2,p);

        % pick the one with the minimum difference
        if d1 < d2
            B = [B(1:I(1)-1,:) ; p ; B(I(1):end,:)];
        else
            B = [B(1:I(1),:) ; p ; B(I(1)+1:end,:)];
        end
    end
    
    if strcmpi(UserData.shape,'bspline')
        B = [B ; B(1,:)];
    end
else
    if I(1) == 1
        % the first point, inject a new point inbetween 1+1 and 1+2
        B = [B(1,:) ; p ; B(2:end,:)];
    elseif I(1) == N
        % the last point, inject a new point inbetween N and N+1
        B = [B(1:N-1,:) ; p ; B(N,:)];
    else
        % was the mouse closer to the next point or the previous point
        Ip = [I(1)-1, I(1)+1];
        
        % compare the distance to the line
        L1 = B([I(1) Ip(1)],:);
        d1 = linepointdistance(L1,p);
        L2 = B([I(1) Ip(2)],:);
        d2 = linepointdistance(L2,p);
        
        % pick the one with the minimum difference
        if (d1 < d2) 
            B = [B(1:I(1)-1,:) ; p ; B(I(1):end,:)];
        else
            B = [B(1:I(1),:) ; p ; B(I(1)+1:end,:)];
        end
    end
end

if strcmpi(UserData.shape,'bspline')
    p = UserData.p;
    Nb = size(B,1);
    t = linspace(0,1,10*Nb-1);
    if UserData.Closed
        % Closed curve
        knots = linspace(0,1,Nb);
        dk = diff(knots(1:2));
        knots = [-(fliplr(dk:dk:p*dk)), knots, 1+(dk:dk:p*dk)];
        phi = bsplines(t,knots,p);
        for k = 1:p
            phi(:,k) = phi(:,k)+phi(:,end-(p-k));
        end
        phi = phi(:,1:Nb-1);
        B = B(1:end-1,:);
    else
        % Open Curve
        Nk = Nb - p + 1;
        knots = [zeros(1,p), linspace(0,1,Nk), ones(1,p)] ;
        phi = bsplines(t,knots,p);
    end
    UserData.phi = phi;
end
B = [mean(B) ; B];

% update the curve coordinates
[x, y, B] = UserData.shapefun(B,UserData);

% update the userdata for possible changes in B
UserData.B = B;
set(H,'UserData',UserData);

% change the coordinates of the plotted curves
set(UserData.h1,'XData',x,'YData',y);
set(UserData.h2,'XData',x,'YData',y);
if UserData.Closed
    set(UserData.h3,'XData',B([2:end, 2],1),'YData',B([2:end, 2],2));
else
    set(UserData.h3,'XData',B(2:end,1),'YData',B(2:end,2));
end
set(UserData.h4,'XData',B(1,1),'YData',B(1,2));

% =====================================================
function d = linepointdistance(L,P)
L2 = (L(2,1)-L(1,1))^2 + (L(2,2)-L(1,2))^2;
if L2 == 0
    d = hypot(L(1,1)-P(1),L(1,2)-P(2));
    return
end
t = ((P(1) - L(1,1)) * (L(2,1) - L(1,1)) + (P(2) - L(1,2)) * (L(2,2) - L(1,2))) / L2;
if (t < 0)
    d = hypot(L(1,1)-P(1),L(1,2)-P(2));
elseif (t > 1) 
    d = hypot(L(2,1)-P(1),L(2,2)-P(2));
else
    d = abs( (L(2,2)-L(1,2))*P(1) - (L(2,1)-L(1,1))*P(2) + L(2,1)*L(1,2) - L(1,1)*L(2,2));
    d = d / sqrt( (L(2,1)-L(1,1))^2 + (L(2,2)-L(1,2))^2 );
end

% =====================================================
function delcontrolpoint(hObject, varargin)
H = varargin{2};
UserData = get(H,'UserData');

% get the current location
p = get(UserData.ha,'CurrentPoint');
p = p(1,1:2);

B = UserData.B;
R = UserData.Rsel;

B = B(2:end,:);

% distance from the mouse to each control point
D = hypot(B(:,1)-p(1),B(:,2)-p(2));
if ~any(D <= R)
    return
end
[D, I] = min(abs(D-R));
N = length(I);

if UserData.Closed == true
    if any(I == [1, N])
        B = B(2:end,:);
        B(end,:) = B(1,:);
    else
        B(I,:) = [];
    end
    if strcmpi(UserData.shape,'bspline')
        B = [B ; B(1,:)];
    end
else
    B(I,:) = [];
end

if strcmpi(UserData.shape,'bspline')
    p = UserData.p;
    Nb = size(B,1);
    t = linspace(0,1,10*Nb-1);
    if UserData.Closed
        % Closed curve
        knots = linspace(0,1,Nb);
        dk = diff(knots(1:2));
        knots = [-(fliplr(dk:dk:p*dk)), knots, 1+(dk:dk:p*dk)];
        phi = bsplines(t,knots,p);
        for k = 1:p
            phi(:,k) = phi(:,k)+phi(:,end-(p-k));
        end
        phi = phi(:,1:Nb-1);
        B = B(1:end-1,:);
    else
        % Open Curve
        Nk = Nb - p + 1;
        knots = [zeros(1,p), linspace(0,1,Nk), ones(1,p)] ;
        phi = bsplines(t,knots,p);
    end
    UserData.phi = phi;
end
B = [mean(B) ; B];
 
% update the curve coordinates
[x, y, B] = UserData.shapefun(B,UserData);

% update the userdata for possible changes in B
UserData.B = B;
set(H,'UserData',UserData);

% change the coordinates of the plotted curves
set(UserData.h1,'XData',x,'YData',y);
set(UserData.h2,'XData',x,'YData',y);
if UserData.Closed
    set(UserData.h3,'XData',B([2:end, 2],1),'YData',B([2:end, 2],2));
else
    set(UserData.h3,'XData',B(2:end,1),'YData',B(2:end,2));
end
set(UserData.h4,'XData',B(1,1),'YData',B(1,2));


% =====================================================
function freemove(hObject, varargin)
% this function controls what happens when nothing is selected. This
% function merely indicates when the mouse is close enough
UserData = get(hObject,'UserData');
p = get(UserData.ha,'CurrentPoint');
p = p(1,1:2);

B = UserData.B;
R = UserData.Rsel;

% distance from the mouse to each control point
D = hypot(B(:,1)-p(1),B(:,2)-p(2));
if any(D <= R)
    % the mouse is hovering a control point
    set(UserData.hf,'Pointer','fleur');
else
    % the mouse is not hovering a control point
    set(UserData.hf,'Pointer','crosshair');
end


% =====================================================
function movehandle(hObject, varargin)
% this function is used after clicking on a control point. It allows the
% modification of a shape.
UserData = get(hObject,'UserData');
if isempty(UserData.CurrentHandle)
    return
end
p = get(UserData.ha,'CurrentPoint');
p = p(1,1:2);

B = UserData.B;
I = UserData.CurrentHandle;

if I == 1;
    % the first control point is always rigid body motion, i.e. move all
    % control points with the same motion
    U = p - B(1,:);
    B(:,1) = B(:,1) + U(1);
    B(:,2) = B(:,2) + U(2);
else
    % any other control point is moved individually
    B(I,:) = p;
end

% update the curve coordinates
[x, y, B] = UserData.shapefun(B,UserData);

% update the userdata for possible changes in B
UserData.B = B;
set(hObject,'UserData',UserData);

% change the coordinates of the plotted curves
set(UserData.h1,'XData',x,'YData',y);
set(UserData.h2,'XData',x,'YData',y);
if UserData.Closed
    set(UserData.h3,'XData',B([2:end, 2],1),'YData',B([2:end, 2],2));
else
    set(UserData.h3,'XData',B(2:end,1),'YData',B(2:end,2));
end
set(UserData.h4,'XData',B(1,1),'YData',B(1,2));

% =====================================================
function selecthandle(hObject, varargin)
% what to do on a mouse click when no control point is currently selected
% (but the mouse can be hovering one)
UserData = get(hObject,'UserData');
m_type = get(hObject, 'selectionType');
if strcmp(m_type, 'open')
    % Detect a double click, exit the tool
    confirmselect([],[],hObject);
end

p = get(UserData.ha,'CurrentPoint');
p = p(1,1:2);

B = UserData.B;
R = UserData.Rsel;

% distance from the mouse to each control point
D = hypot(B(:,1)-p(1),B(:,2)-p(2));
if any(D <= R)
    % the mouse is hovering a control point
    
    % find out which control point
    I = find(D <= R,1);
    if numel(I) == 1
        UserData.CurrentHandle = I;
    else
        [D, I] = min(abs(D-R));
        UserData.CurrentHandle = I;
    end
else
    % the mouse is not hovering a control point
    UserData.CurrentHandle = [];
end

% update the userdata
set(hObject,'UserData',UserData);

if isempty(UserData.CurrentHandle)
    return
end
% change to control point moving mode
set(hObject,'WindowButtonMotionFcn',@movehandle)

% =====================================================
function confirmhandle(hObject, varargin)
% what to do on a mouse click when a control point is currently selected
m_type = get(hObject, 'selectionType');
if strcmp(m_type, 'open')
    % Detect a double click, exit the tool
    confirmselect([],[],hObject);
end
% change to handle free mode
set(hObject,'WindowButtonMotionFcn',@freemove)

% =====================================================
function changecolor(varargin)
% change the line color of the curve
h1 = varargin{3};
h2 = varargin{4};
h3 = varargin{5};
h4 = varargin{6};
col = varargin{7};
set(h1,'Color',col{1});
set(h3,'Color',col{1});
set(h4,'Color',col{1});
if any(col{1} == 'wygc')
    % for light colors use a black background
    set(h2,'Color','k');
    set(h3,'MarkerFaceColor','k');
    set(h4,'MarkerFaceColor','k');
else
    % for dark colors use a white background
    set(h2,'Color','w');
    set(h3,'MarkerFaceColor','w');
    set(h4,'MarkerFaceColor','w');
end

% =====================================================
function confirmselect(varargin)
% stop the interactive mode (i.e. delete the plot h4).
hf = varargin{3};
UserData = get(hf,'UserData');
delete(UserData.h1)
delete(UserData.h2)
delete(UserData.h3)
delete(UserData.h4)


% =====================================================
function keyPressFcn(varargin)
H = varargin{1};
evnt = varargin{2};
UserData = get(H,'UserData');

if isempty(evnt.Modifier)
    evnt.Modifier{1} = '';
end
% fprintf('Key: (%s) %s\n',[evnt.Modifier{:}],evnt.Key);

stride = 0.05;
if any(strcmp(evnt.Modifier,'control'))
    stride = 0.01;
elseif strcmp(evnt.Modifier{1},'shift')
    stride = 0.1;
end

dx = stride*diff(UserData.lim(1,:));
dy = stride*diff(UserData.lim(2,:));
B = UserData.B;


U = [0, 0];
if strcmpi(evnt.Key,'uparrow')
    U = [0, dy];
elseif strcmpi(evnt.Key,'downarrow')
    U = [0, -dy];
elseif strcmpi(evnt.Key,'leftarrow')
    U = [-dx, 0];
elseif strcmpi(evnt.Key,'rightarrow')
    U = [dx, 0];
end
if all(U == 0)
    return
end

if strcmpi(get(UserData.ha,'ydir'),'reverse')
    U(2) = -U(2);
end
if strcmpi(get(UserData.ha,'xdir'),'reverse')
    U(1) = -U(1);
end
    
% apply displacement
B(:,1) = B(:,1) + U(1);
B(:,2) = B(:,2) + U(2);

% update the curve coordinates
[x, y, B] = UserData.shapefun(B,UserData);

% update the userdata for possible changes in B
UserData.B = B;
set(H,'UserData',UserData);

% change the coordinates of the plotted curves
set(UserData.h1,'XData',x,'YData',y);
set(UserData.h2,'XData',x,'YData',y);
if UserData.Closed
    set(UserData.h3,'XData',B([2:end, 2],1),'YData',B([2:end, 2],2));
else
    set(UserData.h3,'XData',B(2:end,1),'YData',B(2:end,2));
end
set(UserData.h4,'XData',B(1,1),'YData',B(1,2));


% =====================================================
function phi = bsplines(x,knots,Np)
% compute bspline shape functions for a given coordinate vector x, knot
% vector knots and bspline order p. Based on the Cox-de Boor recursion
% formula.

n = length(x);
Ni = length(knots);
x = x(:);

% starting with degree zero
p = 0;

% number of basis-functions
Nn = Ni-(p+1);

% initiate matrix for basis functions
phi = zeros(n,Nn);

% Loop over the knots. This defines which coordinates belong to which knot
% section. It is important to use each coordinate only once, and to include
% the edges of the domain in the correct knot sections.
for i = 1:Nn
    if i == (Nn-Np)
        % for the knot on the right edge of x, include the point on the
        % right knot
        I = find(x >= knots(i) & x <= knots(i+1));
    elseif i > (Nn-Np)
        % for the next knots, exclude the point on left knot, include the
        % point on the right knot
        I = find(x > knots(i) & x <= knots(i+1));
    else
        % for all other knots, include the point on the left knot but not
        % on the right knot
        I = find(x >= knots(i) & x < knots(i+1));
    end
    
    % set the function to zero
    phi(I,i) = 1;
end


% Subsequent orders
% =============
for p = 1:Np
    % calculate the number of shape functions for this degree
    Nn = Ni-(p+1);
    
    % store the previous phi as phi1
    phi1 = phi;
    % initiate the current phi (every point supports exactly p+1 shapes if
    % no multiplicity knots are defined, so this is a worst case scenario)
    phi = zeros(n,Nn);
    
    % forloop over the knots
    for i = 1:Nn
        if knots(i+p) == knots(i)
            % if first term == zero (double knot on right side)
            phi(:,i) = phi1(:,i+1).* (knots(i+p+1) - x ) ./ ( knots(i+p+1) - knots(i+1) );
        elseif knots(i+p+1) == knots(i+1)
            % if second term is zero (double knot on left side)
            phi(:,i) = phi1(:,i)  .* (   x - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) ;
        else
            % for all other knots
            phi(:,i) = phi1(:,i)  .* (   x - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) + ...
                phi1(:,i+1).* (knots(i+p+1) - x ) ./ ( knots(i+p+1) - knots(i+1) );
        end
    end
end

phi = sparse(phi);
