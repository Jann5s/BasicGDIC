function [node,conn] = Q8_meshgen(Nrows,Ncols,boundratio,roi)

% Create the nodes
xnode = boundratiovec(roi(1),roi(2),Ncols,boundratio);
ynode = boundratiovec(roi(3),roi(4),Nrows,boundratio);

% place element corner nodes in a larger set
x = zeros(2*Ncols-1,1);
x(1:2:end) = xnode;
y = zeros(2*Nrows-1,1);
y(1:2:end) = ynode;

% place quadratic x-nodes
if Ncols == 2
    x(2) = (x(1)+x(3))/2;
else
    c = xnode(2)-xnode(1);
    d = c/(1 + 1/boundratio);
    d = c/2;
    e = (xnode(3)-xnode(2))/2;
    x(2) = x(1)+d;
    for k = 4:2:length(x)-3
        x(k) = x(k-1)+e;
    end
    x(end-1) = x(end)-d;
end
    
% place quadratic x-nodes
if Nrows == 2
    y(2) = (y(1)+y(3))/2;
else
    c = ynode(2)-ynode(1);
    d = c/(1 + 1/boundratio);
    d = c/2;
    e = (ynode(3)-ynode(2))/2;
    y(2) = y(1)+d;
    for k = 4:2:length(y)-3
        y(k) = y(k-1)+e;
    end
    y(end-1) = y(end)-d;
end

xnode = x;
ynode = y;

% assignin('base','xnode',xnode)
% assignin('base','ynode',ynode)

% if Ncols > 2
%     Lt = (roi(2)-roi(1));
%     Ni = Ncols-2;
%     Li = Lt / (2*boundratio + Ni - 1);
%     Le = boundratio*Li;
%     Lee = Le/(1 + 1/boundratio);
%     xnode = [roi(1) (roi(1)+Lee) (roi(1)+Le):(0.5*Li):(roi(2)-Li) (roi(2)-Lee) roi(2)];
% else
%     xnode = linspace(roi(1),roi(2),3);
% end
% 
% if Nrows > 2
%     Lt = (roi(4)-roi(3));
%     Ni = Nrows-2;
%     Li = Lt / (2*boundratio + Ni - 1);
%     Le = boundratio*Li;
%     Lee = Le/(1 + 1/boundratio);
%     ynode = [roi(3) (roi(3)+Lee) (roi(3)+Le):(0.5*Li):(roi(4)-Li) (roi(4)-Lee) roi(4)];
% else
%     ynode = linspace(roi(3),roi(4),3);
% end

[Xnode, Ynode] = meshgrid(xnode,ynode);
node(:,1) = Xnode(:);
node(:,2) = Ynode(:);

Nnode = size(node,1);
Nel = (Nrows-1)*(Nrows-1);

Inode = reshape(1:Nnode,length(ynode),length(xnode));

conn = zeros(Nel,8);
cnt = 0;
for icol = 1:2:(2*Ncols-2)
    for irow = 1:2:(2*Nrows-2)
        cnt = cnt + 1;
        conn(cnt,:) = [Inode(irow+0,icol  ),... %1
                       Inode(irow+0,icol+2),... %2
                       Inode(irow+2,icol+2),... %3
                       Inode(irow+2,icol  ),... %4
                       Inode(irow+0,icol+1),... %5
                       Inode(irow+1,icol+2),... %6
                       Inode(irow+2,icol+1),... %7
                       Inode(irow+1,icol  )];   %8
    end
end
% used nodes
Iused   = intersect(Inode(:),conn(:));
Iunused = setdiff(Inode(:),conn(:));
for i = 1:length(Iunused)
    A = Iunused(i);
    conn(conn>A) = conn(conn>A) - 1;
    Iunused = Iunused - 1;
end
node = node(Iused,:);

% =======================
