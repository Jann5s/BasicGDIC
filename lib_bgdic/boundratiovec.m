function x = boundratiovec(a,b,N,boundratio)

N = max([N 2]);
N = round(N);
x = zeros(N,1);
boundratio = abs(boundratio);

if N == 2
    x = [a b];
elseif N == 3
    x = [a (a+b)/2 b];
else
    % total length
    Lt = b-a;
    % number of internal nodes
    Ni = N-2;
    % internal element length
    Li = Lt / (2*boundratio + Ni - 1);
    % external element length
    Le = boundratio*Li;
    % knot vector
    x(1) = a;
    x(2) = a+Le;
    for k = 3:(N-1)
        x(k) = x(k-1)+Li;
    end
    x(end) = b;
end
