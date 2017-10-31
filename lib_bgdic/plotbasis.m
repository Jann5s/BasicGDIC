% ==================================================
function [] = plotbasis(varargin)
% plot the current basis
H = varargin{3};
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
% Region of interest
if ~isfield(D,'roi')
    return;
end
roi = D.roi;
% Region of interest
k = get(S.basislist,'Value');
if isfield(D,'basis') && (length(D.basis) >= k) && ~isempty(D.basis(k).name) && isfield(D.basis,'plotphi') && ~isempty(D.basis(k).plotphi)
    basis = D.basis(k);
else
    return;
end

soften = get(S.basissoftenslider,'Value');
soften = round(soften*100)/100;
set(S.basissoftenslider,'Value',soften);
set(S.basissoften,'String',num2str(soften));
alpha = get(S.basisalphaslider,'Value');
alpha = round(alpha*100)/100;
set(S.basisalphaslider,'Value',alpha);
set(S.basisalpha,'String',num2str(alpha));

% get current basis function (i.e. slider position)
id = str2double(get(S.phiid,'String'));
n = D.files(1).size(1);
m = D.files(1).size(2);
x = 1:m;
y = 1:n;
Img = D.files(1).image;
xlim = [1 m];
ylim = [1 n];

% Plot the background image
% =============================================
hi = findobj(H,'Tag','basis_image');
if isempty(hi)
    % color scale
    Img = 255*(Img - min(Img(:))) / (max(Img(:)) - min(Img(:)));
    RGB(:,:,1) = uint8(soften*Img);%+(1-soften)*255;
    RGB(:,:,2) = uint8(soften*Img);%+(1-soften)*255;
    RGB(:,:,3) = uint8(soften*Img);%+(1-soften)*255;
    
    % plot image
    set(0,'CurrentFigure',H)
    set(H,'CurrentAxes',S.axes4)
    set(gca,'NextPlot','replacechildren')
    imagesc(x,y,RGB,'Parent',S.axes4,'Tag','basis_image');
    xlabel('x [px]')
    ylabel('y [px]')
    set(gca,'ydir','reverse')
    set(gca,'xlim',xlim)
    set(gca,'ylim',ylim)
    daspect([1 1 1 ])
else
    Img = 255*(Img - min(Img(:))) / (max(Img(:)) - min(Img(:)));
    RGB(:,:,1) = uint8(soften*Img);%+(1-soften)*255;
    RGB(:,:,2) = uint8(soften*Img);%+(1-soften)*255;
    RGB(:,:,3) = uint8(soften*Img);%+(1-soften)*255;
    set(hi,'CData',RGB);
end

[X, Y] = meshgrid(basis.plotx,basis.ploty);
Iroi = find(X >= roi(1) & X <= roi(2) & Y >= roi(3) & Y <= roi(4));

% plot overlay
% ===============================
plotn = length(basis.ploty);
plotm = length(basis.plotx);
P = reshape(basis.plotphi(:,id),plotn,plotm);
plim(1) = min(P(~isnan(P)));
plim(2) = max(P(~isnan(P)));
if plim(1) == plim(2);
    if plim(1) == 0
        plim = [-0.01, 0.01];
    elseif plim(1) == 1
        plim = [0, 1];
    else
        plim = plim(1) + 0.01*[-plim(1), plim(1)];
    end
end

set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes4)

hi = findobj(H,'Tag','basis_overlay');
if isempty(hi)
    colormap(D.gui.color.cmapoverlay)
    set(gca,'NextPlot','add')
    imagesc(basis.plotx,basis.ploty,P,'Parent',S.axes4,'AlphaData',alpha,'Tag','basis_overlay');
    colorbar;
    hc = colorbar;
    set(hc,'Xcolor',[0 0 0],'Ycolor',[0 0 0],'FontSize',D.gui.fontsize(3))
else
    set(hi,'CData',P,'XData',basis.plotx([1,end]),'YData',basis.ploty([1,end]),'AlphaData',alpha);
end

% corner text
Nphi = basis.Nphi;
str = sprintf('%s, %d/%d',basis.name,id,Nphi);
hi = findobj(H,'Tag','basis_text');
if isempty(hi)
    ht = text(0.02,0.02,str,'units','normalized','Tag','basis_text');
    set(ht,'color',D.gui.color.fg)
    set(ht,'FontSize',D.gui.fontsize(4))
    set(ht,'FontWeight','bold')
else
    set(hi,'String',str);
end
title(str)

% plot mesh
% ===============================
set(gca,'NextPlot','add')
if strcmp(basis.type,'bspline')
    % If mesh type is B-spline
    
    % set the color limits
    clim(1) = 0;
    clim(2) = 1;
    
    [X, Y] = meshgrid(basis.xknot,basis.yknot);
    
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
    
    hi = findobj(H,'Tag','basis_mesh');
    if isempty(hi)
        patch(FV,'FaceColor','none','EdgeColor','k','LineWidth',0.5,'Marker','.','Tag','basis_mesh');
    else
        set(hi,'EdgeColor','k','FaceColor','none','Vertices',FV.Vertices,'Faces',FV.Faces,'Marker','.');
    end
elseif strcmp(basis.type,'zernike')
    % set the color limits
    clim = plim;
    
    % plot the zernike circle
    w = roi(2)-roi(1);
    h = roi(4)-roi(3);
    R = 0.5*sqrt(w^2 + h^2);
    xc = mean(roi(1:2));
    yc = mean(roi(3:4));
    tz = linspace(0,2*pi,60);
    X = R*cos(tz) + xc;
    Y = R*sin(tz) + yc;
    FV.Vertices = [X(:),Y(:)];
    FV.Faces = 1:60;
    
    hi = findobj(H,'Tag','basis_mesh');
    if isempty(hi)
        patch(FV,'FaceColor','none','EdgeColor','k','LineWidth',0.5,'Marker','none','Tag','basis_mesh');
    else
        set(hi,'EdgeColor','k','FaceColor','none','Vertices',FV.Vertices,'Faces',FV.Faces,'Marker','None');
    end
    
    % corner text
    Zp = id - 1;
    Zn = floor(sqrt(Zp)); % Pseudo Zernike
    Zm = Zp - Zn.^2 - Zn; % Pseudo Zernike
    % Zn = ceil((-3+sqrt(9+8*Zp))/2); % Zernike
    % Zm = 2*Zp - Zn.*(Zn+2); % Zernike
    str = sprintf('n = %d, m = %d',Zn,Zm);
    hi = findobj(H,'Tag','basis_text');
    str = [get(hi,'String') ' -- ' str];
    set(hi,'String',str);
    
elseif strcmp(basis.type,'FEM-T')
    % If mesh type is FEM-T
    FV.Vertices = basis.coordinates;
    FV.Faces = basis.connectivity;
    order = basis.order;
    clim = [0,1];
    if order == 2
        clim = [-2/8,1];
        Ipatch = [1 4 2 5 3 6];
        FV.Faces = FV.Faces(:,Ipatch);
    end
    
    hi = findobj(H,'Tag','basis_mesh');
    if isempty(hi)
        patch(FV,'FaceColor','none','EdgeColor','k','LineWidth',0.5,'Marker','s','Tag','basis_mesh');
    else
        set(hi,'EdgeColor','k','FaceColor','none','Vertices',FV.Vertices,'Faces',FV.Faces,'Marker','s');
    end
elseif strcmp(basis.type,'FEM-Q')
    % If mesh type is FEM-Q
    FV.Vertices = basis.coordinates;
    FV.Faces = basis.connectivity;
    order = basis.order;
    clim = [0,1];
    if order == 2
        clim = [-2/8,1];
        Ipatch = [1 5 2 6 3 7 4 8];
        FV.Faces = FV.Faces(:,Ipatch);
    end
    
    hi = findobj(H,'Tag','basis_mesh');
    if isempty(hi)
        patch(FV,'FaceColor','none','EdgeColor','k','LineWidth',0.5,'Marker','s','Tag','basis_mesh');
    else
        set(hi,'EdgeColor','k','FaceColor','none','Vertices',FV.Vertices,'Faces',FV.Faces,'Marker','s');
    end
else
    clim = [-1,1];
    % hide the mesh
    hi = findobj(H,'Tag','basis_mesh');
    if ~isempty(hi)
        set(hi,'EdgeColor','none','Marker','none');
    end
end

% plot the evaluation area
roi = D.roi;
X = [roi(1) roi(2) roi(2) roi(1)];
Y = [roi(3) roi(3) roi(4) roi(4)];
FV.Vertices = [X(:),Y(:)];
FV.Faces = 1:4;

set(0,'CurrentFigure',H)
set(H,'CurrentAxes',S.axes4)
hi = findobj(H,'Tag','basis_roi');
if isempty(hi)
    patch(FV,'FaceColor','none','FaceColor','none','EdgeColor','k','LineWidth',0.5,'Marker','none','Tag','basis_roi');
else
    set(hi,'Vertices',FV.Vertices,'Faces',FV.Faces,'Marker','none');
end

caxis(clim);


drawnow
