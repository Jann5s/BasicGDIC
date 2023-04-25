% ==================================================
function [] = bcwaitbar(H,varargin)
% Button Callback: update the waitbar
% bcwaitbar(H,A)
% H = gui handle
% A = relative part of the waitbar which is full
S = guihandles(H);
D = guidata(H);

if isfield(D,'outputfile')
    % headless mode, don't plot anything
    return
end

C1 = D.gui.color.waitbar1;
C2 = D.gui.color.waitbar2;

% get patch handle
hp = findobj(S.waitbar,'Type','patch');
if ~isempty(hp)
    delete(hp)
end

% get text handle
ht = findobj(S.waitbar,'Type','text');

if nargin == 1
    % only H is specified, clear the wait bar
    if ~isempty(ht)
        delete(ht)
    end
end

if nargin >= 2
    % H and A are specified, update the waitbar to show the status
    A = varargin{1};
    A = max([0 A]);
    A = min([1 A]);
    set(0,'CurrentFigure',H);
    set(H,'CurrentAxes',S.waitbar);
    C(1,1,:) = C1;
    C(1,2,:) = C2;
    X = [0 A;A 1;A 1;0 A];
    Y = [0 0;0 0;1 1;1 1];
    patch(X,Y,C,'edgecolor','none');
    
    % put the text back on top
    uistack(ht,'top');
end

if nargin == 3
    % replace the text
    if ~isempty(ht)
        delete(ht)
    end
    txt = varargin{2};
    ht = text(0.05,0.55,txt,'units','normalized');
    set(ht,'Fontsize',D.gui.fontsize(3))
    set(ht,'VerticalAlignment','Middle')
    set(ht,'HorizontalAlignment','Left')
end
drawnow;
