% ==================================================
function [] = plotres(varargin)
% compute result figure
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if isfield(D,'outputfile')
    % headless mode
    return
end
if ~isfield(D,'files')
    return;
end
if ~isfield(D,'cor')
    return;
end
if ~isfield(D.cor,'converged')
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes8)
    cla
    return;
end
if length(varargin) == 4
    rescompute([],[],H)
end
if ~isfield(D,'res') || isempty(D.res)
    Prompt = 'Select strain definition';
    Title = Prompt;
    SelMode = 'Single';
    List = get(S.resstraindef,'String');
    [selection,ok] = listdlgjn(Prompt,Title,List,SelMode,dlgposition(H));
    if ok == 1
        set(S.resstraindef,'Value',selection)
        % ask for this step
        rescompute([],[],H)
        D = guidata(H);
    else
        return
    end
end
res = D.res;
cor = D.cor;

% options
dimensions = get(S.dicdimensions,'Value')+1; % i.e. 2 or 3 for 2D or 3D

% get the slider positions (and put the rounded values back to the gui)
soften = get(S.ressofteningslider,'Value');
soften = round(soften*100)/100;
set(S.ressofteningslider,'Value',soften);
set(S.ressoftening,'String',num2str(soften));
alpha = get(S.resalphaslider,'Value');
alpha = round(alpha*100)/100;
set(S.resalphaslider,'Value',alpha);
set(S.resalpha,'String',num2str(alpha));
ascale = get(S.resarrowscaleslider,'Value');
ascale = round(ascale*100)/100;
set(S.resarrowscaleslider,'Value',ascale);
set(S.resarrowscale,'String',num2str(ascale));

% get the popupmenu values
str = get(S.resunderlay,'String');
underlay = str{get(S.resunderlay,'Value')};
str = get(S.resoverlay,'String');
overlay = str{get(S.resoverlay,'Value')};
str = get(S.resarrows,'String');
arrows = str{get(S.resarrows,'Value')};

% get the pixel size and unit
pixelsize = str2double(get(S.respixelsize,'String'));
str = get(S.resunit,'String');
unit = str{get(S.resunit,'Value')};

% image size and increment number
n = D.files(1).size(1);
m = D.files(1).size(2);
xlim = [1 m];
ylim = [1 n];
inc = str2double(get(S.resid,'String'));
Ninc = length(cor);
Ncg = 4;
if inc == 0
    inc = 1;
    set(S.resid,'String','1')
    set(S.resslider,'Value',1)
end
if inc > Ninc
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes8)
    cla
    set(gca,'color','none')
    return
end
if isempty(D.res(inc).Exx)
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes8)
    cla
    set(gca,'color','none')
    return
end
if D.cor(inc).done ~= Ncg
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes8)
    cla
    set(gca,'color','none')
    return
end

% plot underlay
% ===============================
if strcmp(underlay,'f')
    A = D.files(1).image;
    A = 255*(A - min(A(:))) / (max(A(:)) - min(A(:)));
    x = 1:m;
    y = 1:n;
elseif strcmp(underlay,'g')
    A = D.files(inc+1).image;
    A = 255*(A - min(A(:))) / (max(A(:)) - min(A(:)));
    x = 1:m;
    y = 1:n;
elseif strcmp(underlay,'gtilde')
    A = cor(inc).froi - cor(inc).r;
    A = 255*(A - min(A(:))) / (max(A(:)) - min(A(:)));
    x = cor(inc).xroi;
    y = cor(inc).yroi;
elseif strcmp(underlay,'r')
    A = cor(inc).r;
    A = 255*(A - min(A(:))) / (max(A(:)) - min(A(:)));
    x = cor(inc).xroi;
    y = cor(inc).yroi;
elseif strcmp(underlay,'none')
    A = ones(n,m)*NaN;
    x = 1:m;
    y = 1:n;
end

% set unit
x = x * pixelsize;
y = y * pixelsize;

set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes8)
xlabel(sprintf('x [%s]',unit))
ylabel(sprintf('y [%s]',unit))
set(gca,'xlim',pixelsize*xlim)
set(gca,'ylim',pixelsize*ylim)
if ~strcmp(underlay,'none')
    % store the underlay as RGB (to fake grayvalue image) to allow using a
    % second colormap
    RGB(:,:,1) = uint8(soften*A);%+(1-soften)*255;
    RGB(:,:,2) = uint8(soften*A);%+(1-soften)*255;
    RGB(:,:,3) = uint8(soften*A);%+(1-soften)*255;
else
    % white background (single pixel rgb image)
    RGB = ones(1,1,3);
end

% plot the RGB
hi = findobj(H,'Tag','res_image');
if isempty(hi)
    set(gca,'NextPlot','replacechildren')
    imagesc(x,y,RGB,'Parent',S.axes8,'Tag','res_image');
    set(gca,'ydir','reverse')
    daspect([1 1 1 ])
else
    set(hi,'CData',RGB,'XData',x([1,end]),'YData',y([1,end]));
end


% plot overlay
% ===============================
plotoverlay = true;
xroi = cor(inc).xroi;
yroi = cor(inc).yroi;
nroi = length(yroi);
mroi = length(xroi);

% set pixel size
xroi = xroi * pixelsize;
yroi = yroi * pixelsize;
[Xroi, Yroi] = meshgrid(xroi,yroi);
Zroi = 0.5*ones(size(Xroi));

U1 = pixelsize*cor(inc).U1;
U2 = pixelsize*cor(inc).U2;
if isfield(cor,'U3') && ~isempty(cor(inc).U3)
    U3 = cor(inc).U3;
else
    U3 = zeros(nroi,mroi);
end
if isfield(cor,'U4') && ~isempty(cor(inc).U4)
    U4 = cor(inc).U4;
else
    U4 = zeros(nroi,mroi);
end
if isfield(cor,'U5') && ~isempty(cor(inc).U5)
    U5 = cor(inc).U5;
else
    U5 = zeros(nroi,mroi);
end
if isfield(cor,'U6') && ~isempty(cor(inc).U6)
    U6 = cor(inc).U6;
else
    U6 = zeros(nroi,mroi);
end

if strcmp(overlay,'U1 (x)')
    B = pixelsize*cor(inc).U1;
    cstr = sprintf('U1 [%s]',unit);
elseif strcmp(overlay,'U2 (y)')
    B = pixelsize*cor(inc).U2;
    cstr = sprintf('U2 [%s]',unit);
elseif strcmp(overlay,'U3 (z or Constant)')
    B = U3;
    if dimensions == 3 % 3D
        cstr = sprintf('U3 [%s]',unit);
    else
        B = B * 100;
        cstr = 'U3 (constant relaxation) [%]';
    end
elseif strcmp(overlay,'U4 (Linear)')
    B = U4 * 100;
    cstr = 'U4 (linear relaxation) [%]';
elseif strcmp(overlay,'U5 (Quadratic)')
    B = U5 * 100;
    cstr = 'U5 (quadratic relaxation) [%]';
elseif strcmp(overlay,'U6 (Cubic)')
    B = U6 * 100;
    cstr = 'U6 (cubic relaxation) [%]';
elseif strcmp(overlay,'strain xx')
    B = res(inc).Exx;
    cstr = '\epsilon_{xx} [-]';
elseif strcmp(overlay,'strain yy')
    B = res(inc).Eyy;
    cstr = '\epsilon_{yy} [-]';
elseif strcmp(overlay,'strain xy')
    B = res(inc).Exy;
    cstr = '\epsilon_{xy} [-]';
elseif strcmp(overlay,'strain yx')
    B = res(inc).Eyx;
    cstr = '\epsilon_{yx} [-]';
elseif strcmp(overlay,'strain maj')
    B = res(inc).Emaj;
    cstr = '\epsilon_{maj.} [-]';
elseif strcmp(overlay,'strain min')
    B = res(inc).Emin;
    cstr = '\epsilon_{min.} [-]';
elseif strcmp(overlay,'strain eq.')
    B = res(inc).Eeq;
    cstr = '\epsilon_{eq.} [-]';
elseif strcmp(overlay,'r')
    B = cor(inc).r;
    B = B*100;
    cstr = 'r (residual) [%]';
elseif strcmp(overlay,'q')
    gt = cor(inc).froi - cor(inc).r;
    B = U3 + gt.*U4 + gt.^2.*U5 + gt.^3.*U6;
    B = B*100;
    cstr = 'q (combined brightness correction) [%]';
elseif strcmp(overlay,'none')
    plotoverlay = true;
    B = ones(nroi,mroi)*NaN;
    cstr = '';
end
if ~exist('B','var') || isempty(B)
    % if for some reason the field is undefined show zeros
    B = ones(nroi,mroi)*NaN;
end
B(cor(inc).Imask) = NaN;
B = double(B);

% assignin('base','B',B)

clim = [-1, 1];
if plotoverlay
    % color limits
    clim1 = str2double(get(S.resclim1,'String'));
    clim2 = str2double(get(S.resclim2,'String'));
    if all([clim1 clim2] == 0);
        maxB = max(max(abs(B(~isnan(B)))));
        if isempty(maxB)
            maxB = 1;
        end
        clim = [-1 1] * maxB;
    else
        clim = [clim1 clim2];
    end
    if clim(1) > clim(2)
        clim = clim([2 1]);
    end
    if clim(1) == clim(2);
        clim = clim(1) + [-0.1 0.1];
    end
    clim = gather(clim);

    % prepare the axes
    set(0,'CurrentFigure',H)
    colormap(D.gui.color.cmapoverlay)
    set(H,'CurrentAxes',S.axes8)
    title(cstr,'FontSize',D.gui.fontsize(4));
    view(2)
    hp = findobj(H,'Tag','res_overlay');
    if isempty(hp)
        set(gca,'NextPlot','add')
        hc = colorbar;
        set(hc,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'FontSize',D.gui.fontsize(3))
        % plot the overlay
        if strcmp(underlay,'g')
            hp = surf(Xroi+U1,Yroi+U2,Zroi,B);
        else
            hp = surf(Xroi,Yroi,Zroi,B);
        end
        set(hp,'FaceAlpha',alpha,'EdgeColor','none','Tag','res_overlay')
    else
        if strcmp(underlay,'g')
            set(hp,'CData',B,'XData',Xroi+U1,'YData',Yroi+U2,'FaceAlpha',alpha);
        else
            set(hp,'CData',B,'XData',Xroi,'YData',Yroi,'FaceAlpha',alpha);
        end
    end
else
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes8)
    hp = findobj(H,'Tag','res_overlay');
    if ~isempty(hp)
        set(hp,'CData',NaN,'XData',x([1,end]),'YData',y([1,end]));
    end
end

ht = findobj(H,'Tag','res_text');
if isempty(ht)
    ht = text(0.02,0.02,sprintf('%d/%d',inc,Ninc),'units','normalized','Tag','res_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(4))
    set(ht,'FontWeight','bold')
else
    set(ht,'String',sprintf('%d/%d',inc,Ninc));
end

% plot arrows
% ===============================
% number of arrows (in x and y)
Nq = 15;

% create a grid
Iqn = round(linspace(2,nroi-1,Nq));
Iqm = round(linspace(2,mroi-1,Nq));

% select only those locations which are not masked
[IQm, IQn] = meshgrid(Iqm,Iqn);
Iq = sub2ind([nroi mroi],IQn(:),IQm(:));
Iq = intersect(Iq,cor(inc).Iunmask);
Xroiq = Xroi(Iq);
Yroiq = Yroi(Iq);

% get the displacements
U1q = pixelsize*cor(inc).U1(Iq);
U2q = pixelsize*cor(inc).U2(Iq);

if strcmp(arrows,'U (x,y)')
    C11 = ascale*U1q;
    C12 = ascale*U2q;
    C21 = nan(length(Iq),1);
    C22 = nan(length(Iq),1);
    Xq = Xroiq;
    Yq = Yroiq;
elseif strcmp(arrows,'strain (xx,yy)')
    C11 = ascale*100*res(inc).Exx(Iq);
    C12 = zeros(length(Iq),1);
    C21 = zeros(length(Iq),1);
    C22 = ascale*100*res(inc).Eyy(Iq);
    if strcmp(underlay,'g')
        Xq = Xroiq+U1q;
        Yq = Yroiq+U2q;
    else
        Xq = Xroiq;
        Yq = Yroiq;
    end
elseif strcmp(arrows,'strain (min,maj)')
    C11 = ascale*100*res(inc).Emin(Iq) .* res(inc).Q11(Iq);
    C12 = ascale*100*res(inc).Emin(Iq) .* res(inc).Q12(Iq);
    C21 = ascale*100*res(inc).Emaj(Iq) .* res(inc).Q21(Iq);
    C22 = ascale*100*res(inc).Emaj(Iq) .* res(inc).Q22(Iq);
    if strcmp(underlay,'g')
        Xq = Xroiq+U1q;
        Yq = Yroiq+U2q;
    else
        Xq = Xroiq;
        Yq = Yroiq;
    end
elseif strcmp(arrows,'none')
    C11 = nan(length(Iq),1);
    C12 = nan(length(Iq),1);
    C21 = nan(length(Iq),1);
    C22 = nan(length(Iq),1);
    Xq = Xroiq;
    Yq = Yroiq;
end

% Force the arrows to be on-top of the overlay
Zq = ones(length(Iq),1);
Uz = zeros(length(Iq),1);

hq1 = findobj(H,'Tag','res_vectors1');
hq2 = findobj(H,'Tag','res_vectors2');
if isempty(hq1)
    hq1 = NaN;
end
if isempty(hq2)
    hq2 = NaN;
end
hq = [hq1,hq2];
if ~strcmp(arrows,'none')
    if any(~ishandle(hq))
        set(gca,'NextPlot','add')
        hq = quiver3(Xq,Yq,Zq,C11,C12,Uz,0,'Tag','res_vectors1');
        set(hq,'Color',D.gui.color.fg);
        set(hq,'LineWidth',0.5);
        hq = quiver3(Xq,Yq,Zq,C21,C22,Uz,0,'Tag','res_vectors2');
        set(hq,'Color',D.gui.color.fg2);
        set(hq,'LineWidth',0.5);
    else
        set(hq,'XData',Xq(:),'YData',Yq(:));
        set(hq(1),'UData',C11(:),'VData',C12(:));
        set(hq(2),'UData',C21(:),'VData',C22(:));
    end
else
    if all(ishandle(hq))
        delete(hq)
    end
end

set(S.secind8,'BackgroundColor',D.gui.color.on)

% set color limits last
caxis(clim)
drawnow

