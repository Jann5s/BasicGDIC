function [node,conn] = meshgen_T3a(Nrows,Ncols,boundratio,roi)

% Create the nodes
xnode = boundratiovec(roi(1),roi(2),Ncols,boundratio);
ynode = boundratiovec(roi(3),roi(4),Nrows,boundratio);

[Xn, Yn] = meshgrid(xnode,ynode);
node = [Xn(:) Yn(:)];

if Ncols > 3
    c = 0.5*diff(xnode(2:3));
    xnode = boundratiovec(xnode(2)-c,xnode(end-1)+c,Ncols-1,1);
else
    xnode = [];
end
if Nrows > 3
    c = 0.5*diff(ynode(2:3));
    ynode = boundratiovec(ynode(2)-c,ynode(end-1)+c,Nrows-1,1);
else
    ynode = [];
end

[Xn, Yn] = meshgrid(xnode,ynode);
node = [node ; [Xn(:) Yn(:)]];

% sort the nodes
sortlist = 1e9*node(:,1) + node(:,2);
[sortlist, I] = sort(sortlist);
node = node(I,:);

conn = delaunay(node(:,1),node(:,2));