% ==================================================
function [] = pateval(varargin)
% Button Callback: Evaluate the pattern
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% load previous pattern evaluations
if isfield(D,'pateval');
    pateval = D.pateval;
else
    pateval = [];
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% update status
D.gui.stat = appendstatus(D.gui.stat,'[2] Pattern evaluation started');

% get current frame
id = str2double(get(S.patid,'String'));

% get the image
A = D.files(id).image;
[n, m] = size(A);
x = 1:m;
y = 1:n;
[X, Y] = meshgrid(x,y);

% set the pattern
if D.gui.activepatview ~= 1
    set(S.patshowpat,'BackgroundColor',D.gui.color.hl);
    set(S.patshowACF,'BackgroundColor',D.gui.color.bg);
    plotpattern(H);
    D.gui.activepatview = 1;
end

% Select an area
set(S.patshowpat,'BackgroundColor',D.gui.color.hl);
set(S.patshowACF,'BackgroundColor',D.gui.color.bg);
set(H,'CurrentAxes',S.axes2)
h = title('position the rectangle, confirm with a doubleclick');
set(h,'color',D.gui.color.axeshl);

% reuse rectangle
if isfield(pateval,'box')
    box = pateval(1).box;
    rect = [box(1),box(3) ; box(2),box(4)];
elseif isfield(D,'roi')
    rect = [D.roi(1),D.roi(3) ; D.roi(2),D.roi(4)];
else
    rect = [0.1*m,0.1*n ; 0.9*m,0.9*n];
end

% load interactive rectangle tool
position = selectarea(rect);
% reset the title
set(h,'color','k');
title('');

box(1) = min(position(:,1));
box(2) = max(position(:,1));
box(3) = min(position(:,2));
box(4) = max(position(:,2));
pateval(1).box = box;

% Crop the image
xlim = [ceil(box(1)) floor(box(2))];
ylim = [ceil(box(3)) floor(box(4))];
A = A(ylim(1):ylim(2),xlim(1):xlim(2));
X = X(ylim(1):ylim(2),xlim(1):xlim(2));
Y = Y(ylim(1):ylim(2),xlim(1):xlim(2));
[n, m] = size(A);

% compute scalar measures
imrms = sqrt(mean((A(:).^2)));
imstd = std(A(:));
imrange = max(A(:)) - min(A(:));
immean = mean(A(:));

acfthreshold = 0.5;%1/exp(1);
[zeta, theta, acf] = correlationlength2d(A,acfthreshold);

% fourier filter
fourierorder = 10;
[four, zeta] = fourier(theta,zeta,fourierorder);

[zetalim(1,1), I] = min(zeta);
zetalim(2,1) = theta(I);
[zetalim(1,2), I] = max(zeta);
zetalim(2,2) = theta(I);

% ZOI auto correlation function (normalized)
% =============================
% number of subsets (Nsub x Nsub)
Nsub = D.gui.ini.patnumfacet.value;

% subset edge indices
In = round(linspace(1,n,Nsub+1));
Im = round(linspace(1,m,Nsub+1));

% subset edge locations
subsetedge.X = X(1,Im);
subsetedge.Y = Y(In,1)';

% evaluate each subset
Z = zeros(60,Nsub*Nsub);
T = zeros(60,Nsub*Nsub);
cnt = 0;
bcwaitbar(H);
for kn = 1:Nsub
    for km = 1:Nsub
        cnt = cnt + 1;
        
        % get subset image
        Asub = A(In(kn):In(kn+1),Im(km):Im(km+1));
        Xsub = X(In(kn):In(kn+1),Im(km):Im(km+1));
        Ysub = Y(In(kn):In(kn+1),Im(km):Im(km+1));

        % get correlation contour
        [Z(:,cnt), T(:,cnt)] = correlationlength2d(Asub,acfthreshold);
        
        % fourier filter
        [four, Z(:,cnt)] = fourier(T(:,cnt),Z(:,cnt),fourierorder);

        

        % store subset data
        subset.xc(cnt) = mean(Xsub(1,:));
        subset.yc(cnt) = mean(Ysub(:,1));
        subset.zeta(cnt) = mean(Z(:,cnt));
        subset.std(cnt) = std(Asub(:));
        subset.mean(cnt) = mean(Asub(:));
    end
    bcwaitbar(H,cnt/(Nsub*Nsub));
end
subset.Z = Z;
subset.T = T;

% store evaluation results to data structure
pateval(id).acf = acf.A;
pateval(id).ux = acf.x;
pateval(id).uy = acf.x;
pateval(id).xlim = xlim;
pateval(id).ylim = ylim;
pateval(id).position = position;
pateval(id).imrms = imrms;
pateval(id).imstd = imstd;
pateval(id).imrange = imrange;
pateval(id).immean = immean;
pateval(id).zeta = zeta;
pateval(id).theta = theta;
pateval(id).zetamean = mean(zeta);
pateval(id).zetalim = zetalim;
pateval(id).subset = subset;
pateval(id).subsetedge = subsetedge;

D.pateval = pateval;

% update info
info(1).str = 'image rms';
info(1).val = imrms;
info(2).str = 'image std';
info(2).val = imstd;
info(3).str = 'image mean';
info(3).val = immean;
info(4).str = 'image range';
info(4).val = imrange;
info(5).str = 'avg. zeta';
info(5).val = mean(zeta);
info(6).str = 'min. zeta';
info(6).val = zetalim(1,1);
info(7).str = 'max. zeta';
info(7).val = zetalim(1,2);
D.gui.info = appendinfo(D.gui.info,info);

% set section indicator on
set(S.secind2,'BackgroundColor',D.gui.color.on)

set(D.gui.waithandles,'Enable','on');drawnow

% update status
D.gui.stat = appendstatus(D.gui.stat,'[2] Pattern evaluation done');
bcwaitbar(H);

% update the application data
guidata(H,D);

plotpattern(H);