% ==================================================
function [] = plotiguess(H)
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

% get current image (i.e. slider position)
Nim = length(D.files);
id = str2double(get(S.iguessid,'String'));
n = D.files(1).size(1);
m = D.files(1).size(2);
x = 1:m;
y = 1:n;
xlim = [1 m];
ylim = [1 n];

% original images
Aref = D.files(1).image;
Acur = D.files(id).image;

% processed images
if isfield(D.files,'imageproc')
    Bref = D.files(1).imageproc;
    Bcur = D.files(id).imageproc;
else
    Bref = Aref;
    Bcur = Acur;
end


for k = 1:4
    if k == 1
        A = Aref;
        ha = S.axes41;
        titlestr = sprintf('image:%d %s',1,'original');
        tag = 'iguess_image1';
    elseif k == 2
        A = Acur;
        ha = S.axes42;
        titlestr = sprintf('image:%d %s',id,'original');
        tag = 'iguess_image2';
    elseif k == 3
        A = Bref;
        ha = S.axes43;
        titlestr = sprintf('image:%d %s',1,'processed');
        tag = 'iguess_image3';
    elseif k == 4
        A = Bcur;
        ha = S.axes44;
        titlestr = sprintf('image:%d %s',id,'processed');
        tag = 'iguess_image4';
    end

    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',ha)
    hi = findobj(H,'Tag',tag);
    title(titlestr)
    if isempty(hi)
        colormap(D.gui.color.cmap)
        set(gca,'NextPlot','replacechildren')
        hi = imagesc(x,y,A,'Parent',ha,'Tag',tag);
        xlabel('x [px]')
        ylabel('y [px]')
        set(gca,'ydir','reverse')
        set(gca,'xlim',xlim)
        set(gca,'ylim',ylim)
        daspect([1 1 1 ])
    else
        set(hi,'CData',A);
    end
    
    clim(1) = min(A(:));
    clim(2) = max(A(:));
    if clim(1) == clim(2);
        clim = clim(1) + [-0.01 0.01];
    end
    caxis(clim);
    
    % If mask is set
    % ===============================
    if ~isempty(D.mask) && any(k == [1 3])
        mask = D.mask;
        
        % maskimage
        Im = maskimage(mask,n,m);
        
        % mask border
        Mb = Im & ~image_erode(Im,1);
        alpha = 0.7;
        
        % alphadata
        AlphaData = zeros(n,m);
        AlphaData(Im)  = alpha;  % the part in the mask is semi-transparent
        AlphaData(~Im) = 1;      % the part outside the mask is transparent
        AlphaData(Mb)  = 0.5;    % the border is opaque
        
        set(hi,'AlphaData',AlphaData);
    else
        set(hi,'AlphaData',1)
    end
    
    ht = findobj(H,'Tag','iguess_text');
    if isempty(ht)
        ht = text(0.01,0.03,sprintf('%d/%d',id,Nim),'units','normalized','Tag','iguess_text');
        set(ht,'color',D.gui.color.fg)
        set(ht,'FontSize',D.gui.fontsize(3))
        set(ht,'FontWeight','bold')
    else
        set(ht,'String',sprintf('%d/%d',id,Nim));
    end
    
end

% If ROI is set
% ===============================
if isfield(D,'roi')
    % plot the evaluation area
    roi = D.roi;
    FV.Vertices = [roi([1 2 2 1])',roi([3 3 4 4])'];
    FV.Faces = [1 2 3 4];
    set(0,'CurrentFigure',H)
    linestyle = '-';
    hp = findobj(H,'Tag','iguess_roi');
    if isempty(hp)
        set(H,'CurrentAxes',S.axes41)
        set(gca,'NextPlot','add')
        hp = patch(FV,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'LineStyle',linestyle,'Tag','iguess_roi');
        set(H,'CurrentAxes',S.axes43)
        set(gca,'NextPlot','add')
        hp = patch(FV,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'LineStyle',linestyle,'Tag','iguess_roi');
    else
        set(hp,'Vertices',FV.Vertices,'LineStyle',linestyle);
    end
end

if isfield(D,'iguess') && ~isempty(D.iguess.x)
    iguess = D.iguess;
    if isfield(iguess,'subset')
        % plot subsets
        hp = findobj(H,'Tag','iguess_subsets');
        if isempty(hp)
            set(0,'CurrentFigure',H)
            set(H,'CurrentAxes',S.axes41)
            set(gca,'NextPlot','add')
            patch(iguess.subset.X,iguess.subset.Y,1,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'Tag','iguess_subsets')
            set(H,'CurrentAxes',S.axes43)
            set(gca,'NextPlot','add')
            patch(iguess.subset.X,iguess.subset.Y,1,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'Tag','iguess_subsets')
        end
    end
    
    % Plot the reference
    href = findobj(H,'Tag','iguess_pointsref');
    set(0,'CurrentFigure',H)
    if ~isempty(href);
        set(href,'XData',iguess.x,'YData',iguess.y);
    else
        set(H,'CurrentAxes',S.axes41)
        set(gca,'NextPlot','add')
        plot(iguess.x,iguess.y,'.','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsref')
        set(H,'CurrentAxes',S.axes42)
        set(gca,'NextPlot','add')
        plot(iguess.x,iguess.y,'.','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsref')
        set(H,'CurrentAxes',S.axes43)
        set(gca,'NextPlot','add')
        plot(iguess.x,iguess.y,'.','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsref')
        set(H,'CurrentAxes',S.axes44)
        set(gca,'NextPlot','add')
        plot(iguess.x,iguess.y,'.','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsref')
        drawnow
    end
    
    % Plot the deformed
    hdef = findobj(H,'Tag','iguess_pointsdef');
    set(0,'CurrentFigure',H)
    if ~isempty(hdef);
        set(hdef,'XData',iguess.x+iguess.ux(:,id),'YData',iguess.y+iguess.uy(:,id));
    else
        set(H,'CurrentAxes',S.axes42)
        set(gca,'NextPlot','add')
        plot(iguess.x+iguess.ux(:,id),iguess.y+iguess.uy(:,id),'o','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsdef')
        set(H,'CurrentAxes',S.axes44)
        set(gca,'NextPlot','add')
        plot(iguess.x+iguess.ux(:,id),iguess.y+iguess.uy(:,id),'o','Color',D.gui.color.fg,'Markersize',1.5*D.gui.markersize,'Tag','iguess_pointsdef')
        drawnow
    end
    
    % set section indicator on
    set(S.secind4,'BackgroundColor',D.gui.color.on)
else
    hs = findobj(H,'Tag','iguess_subsets');
    if ~isempty(hs)
        delete(hs);
    end
    href = findobj(H,'Tag','iguess_pointsref');
    if ~isempty(href)
        delete(href);
    end
    hdef = findobj(H,'Tag','iguess_pointsdef');
    if ~isempty(hdef)
        delete(hdef);
    end
end

drawnow