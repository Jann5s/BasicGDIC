function [nodes,conn] = meshgen_T3b(Nrows,Ncols,boundratio,roi)

% Create the nodes
xnode = boundratiovec(roi(1),roi(2),Ncols,boundratio);
ynode = boundratiovec(roi(3),roi(4),Nrows,boundratio);

[Xn, Yn] = meshgrid(xnode,ynode);
nodes = [Xn(:) Yn(:)];

% sort the nodes
sortlist = 1e9*nodes(:,1) + nodes(:,2);
[sortlist, I] = sort(sortlist);
nodes = nodes(I,:);

conn = delaunay(nodes(:,1),nodes(:,2));