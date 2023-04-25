% ==================================================
function [] = plotcor(H)
% plot the current image to figpanel7
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
n = D.files(1).size(1);
m = D.files(1).size(2);
Nim = length(D.files);
% get current image (i.e. slider position)
inc = str2double(get(S.corid,'String'));
id = inc + 1;
x = 1:m;
y = 1:n;
xlim = x([1 end]);
ylim = y([1 end]);

if isfield(D,'cor') && ...
      (inc <= length(D.cor)) && ...
      (inc > 0) && ...
      isfield(D.cor,'r') && ...
      ~isempty(D.cor(inc).r) && ...
      ~isempty(D.cor(inc).done) && ...
      (D.cor(inc).done > 0)

    cor = D.cor;
        
    % plot image
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes7)
    title(sprintf('residual [%%], inc %d',inc))
    colormap(D.gui.color.cmap)
    hi = findobj(H,'Tag','cor_image');
    if isempty(hi)
        set(gca,'nextplot','replacechildren')
        imagesc(cor(inc).xroi,cor(inc).yroi,cor(inc).r*100,'parent',S.axes7,'Tag','cor_image');
        set(gca,'nextplot','add');
        set(gca,'color','k')
        set(gca,'ydir','reverse')
        set(gca,'xlim',xlim)
        set(gca,'ylim',ylim)
        daspect([1 1 1])
        xlabel('x [px]')
        ylabel('y [px]')
        hc = colorbar;
        set(hc,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'FontSize',D.gui.fontsize(3))
    else
        set(hi,'CData',cor(inc).r*100,'XData',cor(inc).xroi([1,end]),'YData',cor(inc).yroi([1,end]));
        set(gca,'xlim',xlim)
        set(gca,'ylim',ylim)
    end
    
    clim = mean(cor(inc).r(:)) + [-3 3]*std(cor(inc).r(:));
    clim = clim*100;
    caxis(gather(clim));

    % Plot the quiver (arrows)
    Nq = D.gui.ini.cornquiver.value;
    nroi = length(cor(inc).yroi);
    mroi = length(cor(inc).xroi);
    Iqn = round(linspace(2,nroi-1,Nq));
    Iqm = round(linspace(2,mroi-1,Nq));
    [IQm, IQn] = meshgrid(Iqm,Iqn);
    Iq = sub2ind([nroi mroi],IQn(:),IQm(:));
    Iq = intersect(Iq,cor(inc).Iunmask);
    [X, Y] = meshgrid(cor(inc).xroi,cor(inc).yroi);
    hq = findobj(H,'Tag','cor_quiver');
    if isempty(hq)
        hq = quiver(X(Iq),Y(Iq),cor(inc).U1(Iq),cor(inc).U2(Iq),0,'Tag','cor_quiver');
        set(hq,'color',D.gui.color.fg);
    else
        set(hq,'XData',X(Iq),'YData',Y(Iq),'UData',cor(inc).U1(Iq),'VData',cor(inc).U2(Iq));
    end
else
    % plot image
    title(sprintf('image %d',id))
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes7)
    hi = findobj(H,'Tag','cor_image');
    colormap(D.gui.color.cmap)
    if isempty(hi)
        set(gca,'nextplot','replacechildren')
        imagesc(x,y,D.files(id).image,'parent',S.axes7,'Tag','cor_image');
        set(gca,'nextplot','add');
        set(gca,'color','k')
        set(gca,'ydir','reverse')
        set(gca,'ydir','reverse')
        set(gca,'xlim',xlim)
        set(gca,'ylim',ylim)
        daspect([1 1 1])
        xlabel('x [px]')
        ylabel('y [px]')
    else
        set(hi,'CData',D.files(id).image,'XData',x([1,end]),'YData',y([1,end]));
    end
    
    clim(1) = min(D.files(id).image(:));
    clim(2) = max(D.files(id).image(:));
    caxis(clim);
    hc = colorbar;
    set(hc,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'FontSize',D.gui.fontsize(3))
    
    hq = findobj(H,'Tag','cor_quiver');
    if ~isempty(hq)
        X = get(hq,'XData');
        Y = get(hq,'YData');
        U = get(hq,'UData');
        U(:) = NaN;
        set(hq,'XData',X,'YData',Y,'UData',U,'VData',U);
    end
    
end

ht = findobj(H,'Tag','cor_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',id,Nim-1),'units','normalized','Tag','cor_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(3))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',id,Nim-1));
end
drawnow