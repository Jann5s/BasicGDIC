% ==================================================
function [] = plotpattern(H)
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

% hide the ACF contour
delete(findobj(H,'Tag','pat_contour'));

Nim = length(D.files);
% get current frame
id = str2double(get(S.patid,'String'));
n = D.files(1).size(1);
m = D.files(1).size(2);
A = D.files(id).image;

hi = findobj(H,'Tag','pat_image');
if isempty(hi)
    % plot the pattern
    colormap(D.gui.color.cmap)
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes2)
    set(gca,'NextPlot','replacechildren')
    imagesc(A,'Parent',S.axes2,'Tag','pat_image');
    xlabel('x [px]')
    ylabel('y [px]')
    title(sprintf('image:%d name:%s',id,D.files(id).name),'Interpreter','none')
    set(gca,'ydir','reverse')
    set(gca,'xlim',[1,m])
    set(gca,'ylim',[1,n])
    daspect([1 1 1 ])
else
    set(H,'CurrentAxes',S.axes2)
    set(hi,'CData',A,'XData',[1,m],'YData',[1,n]);
    title(sprintf('image:%d name:%s',id,D.files(id).name),'Interpreter','none')
    set(gca,'xlim',[1,m])
    set(gca,'ylim',[1,n])
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

% History Plot
if strcmp(D.files(id).class,'uint8')
    bins = (2^8)-1;
elseif strcmp(D.files(id).class,'uint16')
    bins = (2^16)-1;
else
    bins = (2^8)-1;
end
A = D.files(id).image;
set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes2hist)
cla(S.axes2hist)
set(gca,'NextPlot','replacechildren')
hist(A(:),bins);

% If ROI is set
% ===============================
if isfield(D,'roi')
    % plot the evaluation area
    roi = D.roi;
    FV.Vertices = [roi([1 2 2 1])',roi([3 3 4 4])'];
    FV.Faces = [1 2 3 4];
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes2)
    hp = findobj(H,'Tag','pat_roi');
    if isempty(hp)
        set(gca,'NextPlot','add')
        hp = patch(FV,'EdgeColor',D.gui.color.fg,'FaceColor','none','LineWidth',D.gui.linewidth,'Tag','pat_roi');
    else
        set(hp,'Vertices',FV.Vertices);
    end
end

% If evaluation is done
% ===============================
if ~isfield(D,'pateval')
    return;
end
if length(D.pateval) < id
    return;
end
if isempty(D.pateval(id).acf)
    delete(findobj(H,'Tag','pat_contours'));
    delete(findobj(H,'Tag','pat_boxes'));
    return;
end

% plot subset circles
subset = D.pateval(id).subset;
subsetedge = D.pateval(id).subsetedge;
subsetsize = hypot(mean(diff(D.pateval(id).subsetedge.X)),mean(diff(D.pateval(id).subsetedge.Y)));

% (the data only contains half the contour)
Z = [subset.Z; subset.Z];
T = [subset.T; subset.T+pi];
[Np, Nsub] = size(Z);
X = repmat(subset.xc(:)',Np,1) + Z .* cos(T);
Y = repmat(subset.yc(:)',Np,1) + Z .* sin(T);

% disable ugly contours 
I = (max(abs(subset.Z)) > 0.4*subsetsize);
Y(:,I) = NaN;

set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes2)
set(gca,'NextPlot','Add');
hc = findobj(H,'Tag','pat_contours');
if isempty(hc)
    plot(X,Y,'Color',D.gui.color.fg,'LineWidth',0.5,'Tag','pat_contours')
else
    for k = 1:Nsub
        set(hc(k),'XData',X(:,k),'YData',Y(:,k));
    end
end

% plot subset edges
[X, Y] = meshgrid(subsetedge.X,subsetedge.Y);

% create the connectivity
Nn = numel(X);
[n,m] = size(X);
Ne = (n-1)*(m-1);
In = reshape(1:Nn,n,m);
Ie = reshape(1:Ne,n-1,m-1);
conn = zeros(Ne,4);
for ki = 1:n-1
    for kj = 1:m-1
        ke = Ie(ki,kj);
        con = In(ki:ki+1,kj:kj+1);
        conn(ke,:) = con([1 2 4 3]);
    end
end

FV.Vertices = [X(:),Y(:)];
FV.Faces = conn;

hb = findobj(H,'Tag','pat_boxes');
if isempty(hc)
    patch(FV,'FaceColor','None','EdgeColor',D.gui.color.fg,'LineWidth',0.5,'Tag','pat_boxes');
else
    set(hb,'Vertices',FV.Vertices,'Faces',FV.Faces);
end

% update the title
title('Correlation radius per subset (averaged over all directions)')

clim(1) = min(A(:));
clim(2) = max(A(:));
caxis(clim);

drawnow