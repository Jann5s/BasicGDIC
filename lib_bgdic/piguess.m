function [p, polystr] = piguess(H,D,cg,phi12,phi34)
% compute dof from the predefined initial guess points (auto or manual)

polystr{1} = '';
polystr{2} = '';

% test if there is a initial guess
if ~isfield(D,'iguess') || isempty(D.iguess.x)
    p = zeros(sum(cg.Nphi),1);
    return
end

inc = cg.inc;

x  = D.iguess.x;
y  = D.iguess.y;
ux = D.iguess.ux(:,inc+1);
uy = D.iguess.uy(:,inc+1);

% number of initial guess points
Nu = length(x);

% default poly bases
basis_poly0.type = 'polynomial';
basis_poly0.order = 0;
basis_poly1.type = 'polynomial';
basis_poly1.order = 1;
basis_poly2.type = 'polynomial';
basis_poly2.order = 2;
basis_poly3.type = 'polynomial';
basis_poly3.order = 3;
basis_poly4.type = 'polynomial';
basis_poly4.order = 4;

% fit x-y direction
% ==================================

% interpolate the shape functions on the initial guess locations
phi = zeros(Nu,cg.Nphi(1));
for i = 1:cg.Nphi(1)
    P = reshape(phi12(:,i),cg.n,cg.m);
    P = interp2(cg.x,cg.y,P,x,y,'linear',NaN);
    phi(:,i) = P(:);
end

if rank(phi) == cg.Nphi(1)
    % if there are enough initial guess ponts, perform a direct fit
    p1 = phi \ ux;
    p2 = phi \ uy;
else
    % else, first interpolate with a poly
    
    % poly order depends on number of points
    if Nu <= 3
        basis = basis_poly0;
        polystr{1} = '0';
    elseif Nu <= 6
        basis = basis_poly1;
        polystr{1} = '1st';
    elseif Nu <= 10
        basis = basis_poly2;
        polystr{1} = '2nd';
    elseif Nu <= 15
        basis = basis_poly3;
        polystr{1} = '3rd';
    else
        basis = basis_poly4;
        polystr{1} = '4th';
    end
    
    % create polynomial basis (temporary)
    phitmp = phibuild(cg.x,cg.y,D.roi,basis,H);
    
    phi = zeros(Nu,size(phitmp,2));
    for i = 1:size(phitmp,2)
        P = reshape(phitmp(:,i),cg.n,cg.m);
        P = interp2(cg.x,cg.y,P,x,y,'linear',NaN);
        phi(:,i) = P(:);
    end
    
    % fit poly phi
    p1 = phi \ ux;
    p2 = phi \ uy;
    
    % compute displacement field
    Ux = phitmp * p1;
    Uy = phitmp * p2;
    
    % fit original basis to displacement field
    p1 = phi12 \ Ux(:);
    p2 = phi12 \ Uy(:);
    
    
end
p = [p1 ; p2];

% fit z (brightness) direction
% ==================================
if cg.Ndims >= 3
    f = interp2(cg.x,cg.y,cg.f,x,y);
    g = interp2(cg.x,cg.y,cg.g,x+ux,y+uy);
    uz = g - f;

    % interpolate the shape functions on the initial guess locations
    phi = zeros(Nu,cg.Nphi(3));
    for i = 1:cg.Nphi(3)
        P = reshape(phi34(:,i),cg.n,cg.m);
        P = interp2(cg.x,cg.y,P,x,y,'linear',NaN);
        phi(:,i) = P(:);
    end
    
    if rank(phi) == cg.Nphi(3)
        % if there are enough initial guess ponts, perform a direct fit
        p3 = phi \ uz;
    else
        % else, first interpolate with a poly
        
        % poly order depends on number of points
        if Nu <= 3
            basis = basis_poly0;
            polystr{2} = '0';
        elseif Nu <= 6
            basis = basis_poly1;
            polystr{2} = '1st';
        elseif Nu <= 10
            basis = basis_poly2;
            polystr{2} = '2nd';
        elseif Nu <= 15
            basis = basis_poly3;
            polystr{2} = '3rd';
        else
            basis = basis_poly4;
            polystr{2} = '4th';
        end
        
        % create polynomial basis (temporary)
        phitmp = phibuild(cg.x,cg.y,D.roi,basis,H);
        
        phi = zeros(Nu,size(phitmp,2));
        for i = 1:size(phitmp,2)
            P = reshape(phitmp(:,i),cg.n,cg.m);
            P = interp2(cg.x,cg.y,P,x,y,'linear',NaN);
            phi(:,i) = P(:);
        end
        
        % fit poly phi
        p3 = phi \ uz;
        
        % compute displacement field
        Uz = phitmp * p3;
        
        % fit original basis to displacement field
        p3 = phi34 \ Uz(:);
    end
    p = [p1 ; p2 ; p3];
end

% fit contrast (4) direction
% ==================================
if cg.Ndims >= 4
    p4 = zeros(cg.Nphi(4),1);
    p = [p1 ; p2 ; p3 ; p4];
end

% fit quadratic (5) direction
% ==================================
if cg.Ndims >= 5
    p5 = zeros(cg.Nphi(5),1);
    p = [p1 ; p2 ; p3 ; p4 ; p5];
end

% fit cubic (6) direction
% ==================================
if cg.Ndims == 6
    p6 = zeros(cg.Nphi(6),1);
    p = [p1 ; p2 ; p3 ; p4 ; p5 ; p6];
end

assignin('base','p',p)


