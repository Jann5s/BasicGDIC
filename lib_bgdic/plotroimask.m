% ==================================================
function [] = plotroimask(H)
% plot funtion for the evaluated pattern
S = guihandles(H);
D = guidata(H);
if isfield(D,'outputfile')
    % headless mode
    return
end
if ~isfield(D,'files')
    return;
end

Nim = length(D.files);

id = str2double(get(S.roiid,'String'));
n = D.files(1).size(1);
m = D.files(1).size(2);
A = D.files(id).image;
x = 1:m;
y = 1:n;

xlim = [1 m];
ylim = [1 n];

set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes3)
title('ROI and Mask')

hi = findobj(H,'Tag','roi_image');
if isempty(hi)
    % plot the pattern
    colormap(D.gui.color.cmap)
    set(gca,'NextPlot','replacechildren')
    hi = imagesc(x,y,A,'Parent',S.axes3,'Tag','roi_image');
    xlabel('x [px]')
    ylabel('y [px]')
    set(gca,'ydir','reverse')
    set(gca,'xlim',xlim)
    set(gca,'ylim',ylim)
    daspect([1 1 1 ])
else
    set(hi,'CData',A);
end

ht = findobj(H,'Tag','roi_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',id,Nim),'units','normalized','Tag','roi_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(4))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',id,Nim));
end


% If mask is set
% ===============================
if ~isempty(D.mask)
    mask = D.mask;
    
    % maskimage
    Im = maskimage(mask,n,m);
    
    % mask border
    Mb = Im & ~image_erode(Im,1);
    
    if id == 1
        alpha = 0.6;
    else
        alpha = 0.8;
    end
    
    % alphadata
    AlphaData = zeros(n,m);
    AlphaData(Im)  = alpha;  % the part in the mask is semi-transparent
    AlphaData(~Im) = 1;      % the part outside the mask is transparent
    AlphaData(Mb)  = 0.3;    % the border is opaque
    
    set(hi,'AlphaData',AlphaData);
else
    set(hi,'AlphaData',1);
end

% If ROI is set
% ===============================
if isfield(D,'roi')
    % plot the evaluation area
    roi = D.roi;
    FV.Vertices = [roi([1 2 2 1])',roi([3 3 4 4])'];
    FV.Faces = [1 2 3 4];

    if id == 1
        linestyle = '-';
    else
        linestyle = '-.';
    end
    
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes3)
    hp = findobj(H,'Tag','roi_roi');
    if isempty(hp)
        set(gca,'NextPlot','add')
        hp = patch(FV,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'LineStyle',linestyle,'Tag','roi_roi');
    else
        set(hp,'Vertices',FV.Vertices,'LineStyle',linestyle);
    end
    
    % set section indicator on
    set(S.secind3,'BackgroundColor',D.gui.color.on)
    
end

clim(1) = min(A(:));
clim(2) = max(A(:));
caxis(clim)


drawnow
