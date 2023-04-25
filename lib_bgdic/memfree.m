function [free memsize] = memfree
% output the amount of free memory regardless of OS

if ispc
    [user, sys] = memory;
    free = sys.PhysicalMemory.Available;
    memsize = sys.PhysicalMemory.Total;
elseif ismac
    [r, w] = system('vm_stat');
    w = regexp(w,'\n','split');
    
    i = strncmp('Pages free',w,10);
    fblk = str2double(regexp(w{i}, '[0-9]*', 'match'));
    
    i = strncmp('Pages active',w,12);
    ablk = str2double(regexp(w{i}, '[0-9]*', 'match'));

    i = strncmp('Pages inactive',w,14);
    iblk = str2double(regexp(w{i}, '[0-9]*', 'match'));

    i = strncmp('Pages wired',w,11);
    wblk = str2double(regexp(w{i}, '[0-9]*', 'match'));
    
    free = (fblk+iblk)*4096;
    memsize = (fblk+ablk+iblk+wblk)*4096;
elseif isunix
    [r, w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    memsize = stats(1)*1e3;
    free = (stats(3)+stats(end))*1e3;
end