function node = rotatemesh(node,xc,yc,alpha)
% node = rotatemesh(node,xc,yc,alpha)
% rotate nodes around center [xc,yc] with angle [alpha]

% rotate mesh
beta = atan2(node(:,2)-yc,node(:,1)-xc);
R = sqrt( (node(:,1)-xc).^2 + (node(:,2)-yc).^2 );
ux = R .* cos(alpha*(pi/180)+beta) - (node(:,1)-xc);
uy = R .* sin(alpha*(pi/180)+beta) - (node(:,2)-yc);
node = [node(:,1)+ux,node(:,2)+uy];
