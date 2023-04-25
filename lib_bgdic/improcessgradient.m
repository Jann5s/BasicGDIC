function [f, g] = improcessgradient(cg)

% gradient transform
% ===========
if ~strcmp(cg.improcess,'none')
    [dfdx, dfdy] = gradient(cg.f);
    [dgdx, dgdy] = gradient(cg.g);
end

if strcmp(cg.improcess,'grad|xy|');
    f = sqrt(dfdx.^2 + dfdy.^2);
    g = sqrt(dgdx.^2 + dgdy.^2);
elseif strcmp(cg.improcess,'grad|x|');
    f = abs(dfdx);
    g = abs(dgdx);
elseif strcmp(cg.improcess,'grad|y|');
    f = abs(dfdy);
    g = abs(dgdy);
elseif strcmp(cg.improcess,'grad(xy)');
    f = 0.5*(dfdx + dfdy);
    g = 0.5*(dgdx + dgdy);
elseif strcmp(cg.improcess,'grad(x)');
    f = dfdx;
    g = dgdx;
elseif strcmp(cg.improcess,'grad(y)');
    f = dfdy;
    g = dgdy;
else
    f = cg.f;
    g = cg.g;
end
