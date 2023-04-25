% ==================================================
function [] = plotacf(H)
% plot function for the evaluated autocorrelation function
D = guidata(H);
S = guihandles(H);
id = str2double(get(S.patid,'String'));
if isfield(D,'outputfile')
    % headless mode
    return
end
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);
if ~isfield(D,'pateval')
    return;
end
if length(D.pateval) < id
    return;
end
if isempty(D.pateval(id).acf)
    return;
end

% hide the Pattern stuff
delete(findobj(H,'Tag','pat_contours'));
delete(findobj(H,'Tag','pat_boxes'));
delete(findobj(H,'Tag','pat_roi'));

% get the acf
acf = D.pateval(id).acf;
ux = D.pateval(id).ux;
uy = D.pateval(id).uy;
zeta = D.pateval(id).zeta;
theta = D.pateval(id).theta;

hi = findobj(H,'Tag','pat_image');
if isempty(hi)
    % plot the acf
    set(0,'CurrentFigure',H)
    colormap(D.gui.color.cmap)
    set(H,'CurrentAxes',S.axes2)
    set(gca,'NextPlot','replacechildren')
    hi = imagesc(ux,uy,acf,'Parent',S.axes2,'Tag','pat_image');
    colorbar
    xlabel('ux [px]')
    ylabel('uy [px]')
    set(gca,'ydir','reverse')
    set(gca,'xlim',3*D.pateval(id).zetamean*[-1,1])
    set(gca,'ylim',3*D.pateval(id).zetamean*[-1,1])
    daspect([1 1 1 ])
else
    set(H,'CurrentAxes',S.axes2)
    set(hi,'CData',acf,'XData',ux([1,end]),'YData',uy([1,end]));
    set(gca,'xlim',3*D.pateval(id).zetamean*[-1,1])
    set(gca,'ylim',3*D.pateval(id).zetamean*[-1,1])
    caxis('auto');
end

ht = findobj(H,'Tag','pat_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',id,Nim),'units','normalized','Tag','pat_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(4))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',id,Nim));
end

% plot the contour (the data only contains half the contour)
zeta = [zeta; zeta];
theta = [theta ; theta+pi];
X = zeta.*cos(theta);
Y = zeta.*sin(theta);
hc = findobj(H,'Tag','pat_contour');
if isempty(hc)
    set(gca,'NextPlot','add')
    plot(X,Y,'Color',D.gui.color.fg,'LineWidth',D.gui.linewidth,'Tag','pat_contour')
else
    set(hc,'XData',X,'YData',Y);
end

% update the title
title('auto-correlation function with correlation contour')