% ==================================================
function [] = plotimage(H)
% plot the current image to figpanel1
S = guihandles(H);
D = guidata(H);
if isfield(D,'outputfile')
    % headless mode
    return
end
% test if there are images
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get current image (i.e. slider position)
id = str2double(get(S.imageid,'String'));
n = D.files(1).size(1);
m = D.files(1).size(2);
A = D.files(id).image;
[n, m] = size(A);
xlim = [1 m];
ylim = [1 n];
x = 1:m;
y = 1:n;

% plot image
set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes1)
title(sprintf('image:%d name:%s',id,D.files(id).name),'Interpreter','none')

hi = findobj(H,'Tag','image_image');
if isempty(hi)
    colormap(D.gui.color.cmap)
    set(gca,'NextPlot','replacechildren')
    imagesc(x,y,A,'Parent',S.axes1,'Tag','image_image');
    xlabel('x [px]')
    ylabel('y [px]')
    set(gca,'ydir','reverse')
    set(gca,'xlim',xlim)
    set(gca,'ylim',ylim)
    daspect([1 1 1 ])
else
    set(hi,'CData',A);
end

ht = findobj(H,'Tag','image_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',id,Nim),'units','normalized','Tag','image_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(4))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',id,Nim));
end
drawnow