function phi = phibuild(x,y,roi,basis,H,varargin)
% create the basis functions for a given mesh and image

% debug
% ===============
% clear all; close all
% load('mesh.mat','mesh');
% roi = [51.2000  460.8000   38.4000  345.6000];
% n = 384;
% m = 512;
% x = 1:m;
% y = 1:n;
% 
% mesh.type = 'poly3';
% % mesh.type = 'FEM-T';

precision = 'double';
if nargin == 6
    if strcmpi(varargin{1},'single')
        precision = 'single';
    end
end

if strcmp(basis.type,'polynomial')
    % polynomial basis
    p = basis.order;

    n = 0:p;
    cnt = 0;
    for k = n
        N = 0:k;
        for kk = 1:length(N)
            cnt = cnt + 1;
            philist(cnt,:) = [N(1+(kk-1)) N(end-(kk-1)) 1];
        end
    end

    if 0 == 1
        phi = buildphi_poly(x,y,philist,roi,precision);
    elseif 0 == 1
        phi = buildphi_chebyshev(x,y,philist,roi,precision);
    elseif 1 == 1
        phi = buildphi_legendre(x,y,philist,roi,precision);
    end
    phi = phi.x;

elseif strcmp(basis.type,'harmonic')
    % harmonic basis
    p = basis.order;

    phi = buildphi_harmonic(x,y,p,roi,H,precision);

elseif strcmp(basis.type,'zernike')
    % harmonic basis
    p = basis.order;
    
    if false
        phi = buildphi_zernike(x,y,p,roi,H,precision);
    else
        phi = buildphi_pseudozernike(x,y,p,roi,H,precision);
    end
    
    
elseif strcmp(basis.type,'FEM-T')
    % FEM triangles
    nodes = basis.coordinates;
    conn = basis.connectivity;
    order = basis.order;
    
    if order == 1
        phi = buildphi_femt3(x,y,nodes,conn,H,precision);
        phi = phi.a;
    elseif order == 2
        phi = buildphi_femt6(x,y,nodes,conn,H,precision);
    end
    
elseif strcmp(basis.type,'FEM-Q')
    % FEM triangles
    nodes = basis.coordinates;
    conn = basis.connectivity;
    order = basis.order;
    
    if order == 1
        phi = buildphi_femq4(x,y,nodes,conn,H,precision);
    elseif order == 2
        phi = buildphi_femq8(x,y,nodes,conn,H,precision);
    end
    
elseif strcmp(basis.type,'bspline')    
    % B-Spline mesh
    p = basis.order;
    xknot = basis.xknot;
    yknot = basis.yknot;
    
    phi = buildphi_bspline(x,y,xknot,yknot,p,roi,H,precision);
end