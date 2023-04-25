function headlessstatus(str)

if iscell(str)
    for k = 1:length(str)
        headlessstatus(str{k})
    end
else
    
    
    % add status line
    if strcmp(str,'rule')
        str = '=============================================================';
    elseif strcmp(str,'smallrule')
        str = '-----------------------------------';
    end
    
    prefix = datestr(now,'HH:MM:SS');
    fprintf(1,'%s, %s\n',prefix,str);
end
