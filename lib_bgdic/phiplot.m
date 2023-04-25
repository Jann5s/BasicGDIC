function phiplot(x,y,phi,varargin)
%PHIPLOT plots all the basis functions of the matrices found in phi, i.e.
%    phi.x, phi.y, phi.z
%    
%    phiplot(x,y,phi)
%
%    or
%
%    phiplot(x,y,phi,'unit')
%
%    See also JNDIC, BUILDPHI_POLY, BUILDPHI_LEGENDRE, BUILDPHI_CHEBYSHEV, 
%        BUILDPHI_FEMT3.

if length(x) == 1
    m = x;
    x = 1:m;
else
    m = length(x);
end
if length(y) == 1
    n = y;
    y = 1:n;
else
    n = length(y);
end
if nargin == 4
    unit = [' [' varargin{1} ']'];
else
    % set the unit string
    unit = ' [mm]';
end

% set the color limits
clim = [-1 1];
% font options
fontsize = 16;
fontname = 'times';
% figure size and position
figpos = [100 100 600 400];

% test if phi is not a struct
if ~isstruct(phi)
    % create a struct with only one field
    P.a = phi;
    phi = P;
end

% find the important fields (i.e. remove phi.n and phi.m from the list)
fields = fieldnames(phi);
Nf = length(fields);
dirs = {};
for k = 1:Nf
    if any(strcmpi(fields{k},{'a';'x';'y';'z'}))
        dirs = [dirs ; fields{k}];
    end
end

% plot
% ======================
Nd = length(dirs);
for kk = 1:Nd
    N = size(phi.(dirs{kk}),2);

    for k = 1:N
        % rebuild the basis function as an image
        P = reshape(phi.(dirs{kk})(:,k),n,m);
        
        % build some strings
        figname = sprintf('phi.%s(%d)',dirs{kk},k);
        cbarstr = sprintf('phi.%s(%d)',dirs{kk},k);
        
        % create a new fig window
        figure('Name',figname,'Position',figpos);
        % set the colormap
        colormap(jnmap)
        % plot the image
        imagesc(x,y,P);
        % set the color limits
        caxis(clim)
        % set some axes properties
        set(gca,'FontSize',fontsize)
        set(gca,'FontName',fontname)
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'ydir','normal')
        % set the labels
        xlabel(['x' unit])
        ylabel(['y' unit])
        % plot a colorbar
        hc = colorbar('fontsize',fontsize,'FontName',fontname);
        % label the colorbar
        set(get(hc,'YLabel'),'String',cbarstr,'fontsize',fontsize,'FontName',fontname)
    end
end

