% ==================================================
function dlgpos = dlgposition(H)
% helper function to get the center position for the dialog box, used for
% inputdlgjn

guipos = get(H,'Position');
xc = guipos(1) + 0.3*guipos(3);
yc = guipos(2) + 0.8*guipos(4);

dlgpos = [xc yc];
