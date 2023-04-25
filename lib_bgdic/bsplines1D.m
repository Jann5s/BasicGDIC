function N = bsplines1D(x,knots,p)
n = length(x);
Np = p;
x = x(:);
% dof = Ni + (Np-1);

% 1D B-Spline Shape Functions
% ========================


% Make the B-Spline Open: repeat end knots p times
knots = [knots(1)*ones(Np,1) ; knots(:) ; knots(end)*ones(Np,1)];
Ni = length(knots);

% Zero order
% =============

% degree
p = 0;

% number of basis-functions
Nn = Ni-(p+1);

% initiate matrix for basis functions
N0 = zeros(n,Nn);

for i = 1:Nn
    if i == (Nn-Np)
        % for the knot on the right edge of x, include the point on the
        % right knot
        I = find(x >= knots(i) & x <= knots(i+1));
    elseif i > (Nn-Np)
        % for the next knots, exclude the point left knot, include the
        % point on the right knot
        I = find(x > knots(i) & x <= knots(i+1));
    else
        % for all other knots, include the point on the left knot but not
        % on the right knot
        I = find(x >= knots(i) & x < knots(i+1));
    end    
    % set the function to zero
    N0(I,i) = 1;
end
% copy the zero basis for later use
N = N0;

% Subsequent orders
% =============
for p = 1:Np
    % calculate the number of shape functions for this degree
    Nn = Ni-(p+1);
    
    % store the previous N as N1
    N1 = N;
    % initiate the current N
    N  = zeros(n,Nn);
    
    % forloop over the knots
    for i = 1:Nn
        if knots(i+p) == knots(i)
            % if first term == zero (double knot on right side)
            N(:,i) = N1(:,i+1).* ( knots(i+p+1) - x ) ./ ( knots(i+p+1) - knots(i+1) );
        elseif knots(i+p+1) == knots(i+1)
            % if second term is zero (double knot on left side)
            N(:,i) = N1(:,i)  .* (   x - knots(i)   ) ./ ( knots(i+p)   - knots(i)   ) ;
        else
            % for all other knots
            N(:,i) = N1(:,i)  .* (   x - knots(i)   ) ./ ( knots(i+p)   - knots(i)   ) + ...
                     N1(:,i+1).* ( knots(i+p+1) - x ) ./ ( knots(i+p+1) - knots(i+1) );
        end
    end
    
end    
