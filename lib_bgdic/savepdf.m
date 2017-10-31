function savepdf(filename,varargin)

if nargin == 1
    H = gcf;
elseif nargin == 2
    H = varargin{1};
else
    error('savepdf: incorrect number of inputs');
end

% fix the extention, to be always .png (small case)
filename = [regexprep(filename,'.pdf$','','ignorecase') '.pdf'];

% get the original figure position (and size)
savepos = get(H,'Position');
% set the paper position to 1 inch per 100 pixels
set(H,'PaperUnits','inches','PaperPosition',savepos.*[0 0 1e-2 1e-2])
set(H,'PaperSize',savepos(3:4).*[1e-2 1e-2])

% save pdf
print(H,filename,'-dpdf')
