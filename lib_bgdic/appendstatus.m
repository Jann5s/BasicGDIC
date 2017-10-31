function stat = appendstatus(varargin)

if nargin == 1
    stat = {''};
    str = varargin{1};
elseif nargin == 2
    stat = varargin{1};
    str  = varargin{2};
else
    return
end

% add status line
if strcmp(str,'rule')
    str = '=============================================================';
elseif strcmp(str,'smallrule')
    str = '-----------------------------------';
end

prefix = datestr(now,'HH:MM:SS');
str = [prefix ', ' str];
if (length(stat) == 1) && isempty(stat{1})
    stat = {str};
else
    stat = [str ; stat];
end
