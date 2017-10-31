function [f, g] = improcessblur(cg)

% blur
% ===========
if cg.blur > 0.3
    f = image_blur(cg.f,cg.blur);
    g = image_blur(cg.g,cg.blur);
else
    f = cg.f;
    g = cg.g;
end
