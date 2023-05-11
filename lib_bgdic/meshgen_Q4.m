function [node,conn] = Q4_meshgen(Nrows,Ncols,boundratio,roi)

% Create the nodes
xnode = boundratiovec(roi(1),roi(2),Ncols,boundratio);
ynode = boundratiovec(roi(3),roi(4),Nrows,boundratio);

[Xnode, Ynode] = meshgrid(xnode,ynode);
node(:,1) = Xnode(:);
node(:,2) = Ynode(:);

Nnode = size(node,1);
Nel = (Nrows-1)*(Ncols-1);

% Connect the nodes
Inode = reshape(1:Nnode,length(ynode),length(xnode));

conn = zeros(Nel,4);
cnt = 0;
for icol = 1:Ncols-1
    for irow = 1:Nrows-1
        cnt = cnt + 1;
        conn(cnt,:) = [Inode(irow  ,icol  ),...
                       Inode(irow  ,icol+1),...
                       Inode(irow+1,icol+1),...
                       Inode(irow+1,icol  )];
    end
end
