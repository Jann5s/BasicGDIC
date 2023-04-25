function phi = dofbuild_poly(n)
% this functions builds a phi matrix which has all shape functions
% "activated" up to the order n for each direction. The resulting phi
% matrix can be directly supplied to globalDIC2D.m. This function works
% well for the "polynomial" and the "Chebyshev" families of functions.

if length(n)==1
    n = 0:n;
end

cnt = 0;
for k = n
    N = 0:k;
    for kk = 1:length(N)
        cnt = cnt + 1;
        phi(cnt,:) = [N(1+(kk-1)) N(end-(kk-1)) 1];
        cnt = cnt + 1;
        phi(cnt,:) = [N(1+(kk-1)) N(end-(kk-1)) 2];
        cnt = cnt + 1;
        phi(cnt,:) = [N(1+(kk-1)) N(end-(kk-1)) 3];
    end
end
