function P = harmonic1D(x,p)
if p == 0
    dof = 1;
else
    dof = 1 + 2*p;
end

n = length(x);
P = zeros(n,dof);
period = 4;

if true
    % new style
    N = @(x,k) 0.5*(1 + (-1)^k)*sin(k*pi*x/period) + 0.5*(1 + (-1)^(k-1))*cos((k-1)*pi*x/period);
    
    for k = 1:dof
        P(:,k) = N(x,k);
    end
    
else
    % old style
    
    % first shape is just constant
    for k = 1:2:dof
        if k == 1
            P(:,k) = 1;
        else
            P(:,k-1) = sin((k-1)*pi*x/period);
            P(:,k)   = cos((k-1)*pi*x/period);
        end
    end
end
