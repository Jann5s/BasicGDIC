function [node,conn] = meshgen_T3toT6(node,conn,boundratio,roi)
Nel = size(conn,1);

% Add quadratic nodes
% =========================

conn2 = zeros(Nel,6);
conn2(:,1:3) = conn;
          
for iel = 1:Nel
    % node coordinates
    xel = [node(conn(iel,:),1), node(conn(iel,:),2)];
    
    I1 = [1 2 3];
    I2 = [2 3 1];

    for i = 1:3
        % edge nodes
        x1 = xel(I1(i),1);
        x2 = xel(I2(i),1);
        y1 = xel(I1(i),2);
        y2 = xel(I2(i),2);
        
        % edge angle
        alpha = atan2(y2-y1,x2-x1);
        
        % edge length
        R = sqrt( (x2-x1)^2 + (y2-y1)^2 );
        
        Nnode = size(node,1);
        
        % extra node position (depending if an element edge is on the outer
        % boundary of the domain
        if any(x1 == roi(1:2)) && any(y1 == roi(3:4)) && ...
           any(x2 == roi(1:2)) && any(y2 == roi(3:4))
            % first and second node are a corner nodes
            a = 0.5;
        elseif any(x1 == roi(1:2)) && any(y1 == roi(3:4))
            % first node is a corner node
            a = 1/(1+1/boundratio);
        elseif any(x2 == roi(1:2)) && any(y2 == roi(3:4))
            % second node is a corner node
            a = 1/(1+boundratio);
        elseif any(x1 == roi(1:2)) && all(x2 ~= roi(1:2)) && any(y2 == roi(3:4)) && all(y1 ~= roi(3:4))
            % each node is on a different edge but not a corner node (type 1)
            a = 0.5;
        elseif any(x2 == roi(1:2)) && all(x1 ~= roi(1:2)) && any(y1 == roi(3:4)) && all(y2 ~= roi(3:4))
            % each node is on a different edge but not a corner node (type 2)
            a = 0.5;
        elseif (x1 == roi(1)) && (x2 == roi(1))
            % both nodes on left edge
            a = 0.5;
        elseif (x1 == roi(2)) && (x2 == roi(2))
            % both nodes on right edge
            a = 0.5;
        elseif (y1 == roi(3)) && (y2 == roi(3))
            % both nodes on top edge
            a = 0.5;
        elseif (y1 == roi(4)) && (y2 == roi(4))
            % both nodes on bottom edge
            a = 0.5;
        elseif any(x1 == roi(1:2)) || any(y1 == roi(3:4))
            % first node on edge
            a = 1/(1+1/boundratio);
        elseif any(x2 == roi(1:2)) || any(y2 == roi(3:4))
            % second node on edge
            a = 1/(1+boundratio);
        else
            % no nodes on edge
            a = 0.5;
        end
        
        % extra nodes
        node(Nnode+1,1) = a*R * cos(alpha) + x1;
        node(Nnode+1,2) = a*R * sin(alpha) + y1;
        conn2(iel,3+i) = Nnode+1;
    end
end
Nnode = size(node,1);
conn = conn2;

% Sweep (find double nodes)
sweeptol = 1e-6;
Iuniq = 1:Nnode;
for i = 1:Nnode
    I = find( abs(node(Iuniq,1) - node(i,1)) < sweeptol &...
              abs(node(Iuniq,2) - node(i,2)) < sweeptol);
    if length(I) > 1
        for k = 2:length(I)
            conn(conn == I(k)) = I(1);
            node(I(k),:) = NaN;
        end
        
    end
    
    
end

% Remove unused nodes
Inode   = 1:Nnode;
Iused   = intersect(Inode(:),conn(:));
Iunused = setdiff(Inode(:),conn(:));
for i = 1:length(Iunused)
    A = Iunused(i);
    conn(conn>A) = conn(conn>A) - 1;
    Iunused = Iunused - 1;
end
node = node(Iused,:);

% sort the nodes
sortlist = 1e9*node(:,1) + node(:,2);
[sortlist, I] = sort(sortlist);
node = node(I,:);

conn2 = conn;
Nnode = size(node,1);
for k = 1:Nnode
    conn(conn2==I(k)) = k;
end




