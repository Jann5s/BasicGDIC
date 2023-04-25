function msgdlgjn(Message,Position)
% function similar to msgbox, except simpler and with the option to
% position it.

if ~iscell(Message)
    Message = {Message};
end

Title = 'basic gdic message';

fontsize = 12;
N = length(Message);
M = max(cellfun('length',Message));

Figheight = N*2.5*fontsize + 3.0*fontsize;
Figwidth = M*0.8*fontsize + 2*fontsize;
Figpos(3:4) = [Figwidth Figheight];
Figpos(1:2) = Position - 0.5*Figpos(3:4);

FigColor=get(0,'DefaultUicontrolBackgroundColor');
WindowStyle='modal';
Interpreter='none';
Resize = 'off';

H =    dialog(                     ...
  'Visible'          ,'on'      , ...
  'Name'             ,Title      , ...
  'Pointer'          ,'arrow'    , ...
  'Units'            ,'pixels'   , ...
  'UserData'         ,'Cancel'   , ...
  'Tag'              ,Title      , ...
  'HandleVisibility' ,'callback' , ...
  'Color'            ,FigColor   , ...
  'NextPlot'         ,'add'      , ...
  'WindowStyle'      ,WindowStyle, ...
  'Resize'           ,Resize,      ...
  'Position'         ,Figpos       ...
  );


xpos = linspace(0.08,0.92,4);
ypos(1) = 0.08;
ypos(2) = ypos(1) + 0.05 + 2/(2*N+2);
height = 1.8/(1.5*N+2);
width = 0.9*mean(diff(xpos));

% Create title
uicontrol('String',Message,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'units','normalized',...
    'Position',[xpos(1) ypos(2) xpos(end)-xpos(1) 1-(ypos(2)+0.06)],...
    'FontSize',fontsize,...
    'Parent',H);

uicontrol('String','Ok',...
    'ToolTipString','Confirm',...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[xpos(3) ypos(1) width height],...
    'FontSize',fontsize,...
    'Parent',H,...
    'call',{@pressOk,H});


set(H,'KeyPressFcn',{@doFigureKeyPress,H})

uiwait(H)

if ishandle(H)
    delete(H);
end

function doFigureKeyPress(obj, evd, H) %#ok
switch(evd.Key)
  case {'return','space'}
      pressOk([],[],H);
  case {'escape'}
      pressOk([],[],H);
end

function pressOk(varargin)
H = varargin{3};
uiresume(H);


