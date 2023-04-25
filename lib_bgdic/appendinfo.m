function info = appendinfo(iold,inew)
% update info (replace if entry exists)

info = iold;
for k = 1:length(inew)
    str = inew(k).str;
    I = strcmp({info(:).str},str);
    if any(I)
        info(I) = inew(k);
    else
        info(end+1) = inew(k);
    end
end