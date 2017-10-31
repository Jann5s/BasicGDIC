function varargout = basicgdic(varargin)
%   BASICGDIC() is a program with a graphical user interface which can be
%   used to perform Global Digital Image Correlation. Please see the help
%   which is embedded in the user interface in the <Info> tab.
%
%   Run this file to start the program, the help is also contained in the
%   user interface.
%
%   Optional:
%     - basicgdic(inputfile), where inputfile is a string specifying a GUI
%       state file as saved in section 9. This file will be loaded after
%       the GUI starts.
%     - basicgdic(D), where D is a structure as produced by the 
%       <Guidata to D> button found in section 9.
%     - basicgdic(inputfile,savefile), which initates the GUI in an
%       invisible state (for headless computers). The input file will be
%       loaded, the correlation started, and the state file will be saved
%       in the file specified by the savefile string.
%
%   Version: 1.0
%
%   Copyright 2017 Jan Neggers 
%   This program is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%   General Public License for more details.
%
%   Both the General Public License version 3 (GPLv3), and the Lesser
%   General Public License (LGPLv3) are supplied with this program and can
%   be found in the lib_bgdic folder as gpl.txt and lgpl.txt respectively.

bgdicversion = 1.0;
D.version = bgdicversion;

% input testing
headlessmode = false;
if nargin == 1
    inputfile = varargin{1};
elseif nargin == 2
    headlessmode = true;
    inputfile = varargin{1};
    outputfile = varargin{2};
end

% Add library
addpath(fullfile(pwd,'lib_bgdic'));

% get the display size in pixels
screensize = get(0,'screensize');

% GUI size
pos.h = 700  ;    % Fig height
pos.w = 1200 ;    % Fig width
pos.x = (screensize(3)-pos.w)/2;  % Fig left
pos.y = (screensize(4)-pos.h)/2;  % Fig bottom

% define fontsizes
% --------------
fontsize(1) = 9;
fontsize(2) = 11;
fontsize(3) = 13;
fontsize(4) = 18;
if ispc
    fontsize = fontsize - 2;
    fixwidthfont = 'FixedWidth';
else
    fontsize = fontsize - 1;
    fixwidthfont = 'Monospaced';
end

% define colors
% --------------
color.fg     = [1.0, 0.4, 0.0];
color.fg2    = [0.8, 0.0, 0.5];
color.bg     = [0.8, 0.8, 0.8];
color.hl     = [0.6, 0.7, 1.0];
color.on     = [0.0, 0.7, 0.0];
color.init   = [0.7, 0.7, 0.7];
color.off    = [0.9, 0.7, 0.1];
color.axeshl = [0.0, 0.0, 0.7];
color.axesfg = color.bg;
color.waitbar1 = color.hl;
color.waitbar2 = color.fg;
color.c{1,:} = [0.8 0 0  ];
color.c{2,:} = [0.8 0 0  ];
color.c{3,:} = [0   0 0.8];

% colormaps
Nc = 64;
color.cmap = gray(Nc);
color.cmapoverlay = parula(Nc);
% color.cmapoverlay = jet(Nc);

% load the defaults
ini = iniread;

% annoying warnings
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
warning('off','MATLAB:hg:patch:RGBColorDataNotSupported')

% load GUI elements
basicgdic_uicontrols;
drawnow;
set(0,'CurrentFigure',H);

% set the renderer
set(H,'Renderer','opengl')
try
    opengl hardwarebasic
end

% Defaults
% =====================
D.gui.activepanel = 9;
D.gui.activeinfopanel = 3;
D.gui.activepatview = 1;
D.gui.ini = ini;
D.mask = [];
D.basis(4).name = [];

% initialize info
D.gui.info(1).str = 'date';
D.gui.info(1).val = datestr(now,'yyyy-mm-dd');
D.gui.info(2).str = 'time';
D.gui.info(2).val = datestr(now,'HH:MM:SS');

% initialize status
str = ['Basic gdic tool started: ' datestr(now,'yyyy-mm-dd')];
D.gui.stat = appendstatus({''},str);
D.gui.stat = appendstatus(D.gui.stat,'rule');

% initialize help
helpstr = basicgdic_help;
set(Hhelp,'String',helpstr);

% Finalizing GUI buildup
% =====================

% convert to normalized units (to allow resizing)
Hu = findobj(H,'-property','units');
Hu = setdiff(Hu,H);
set(Hu,'units','normalized');

% Generate handles structure
S = guihandles(H);

% add to data structure
D.gui.pos = pos;
D.gui.color = color;
D.gui.fontsize = fontsize;
D.gui.fixwidthfont = fixwidthfont;
D.gui.linewidth = 1.6;
D.gui.markersize = 12;

% collection of controls to disable when processing
D.gui.waithandles = waithandlesfun(S);

% update the application data
guidata(H,D);
basisset(H);
                        
% debugging
% --------------
% assignin('base','H',H)
% assignin('base','S',S)
% assignin('base','D',D)

% read inputfile from ini (the input argument takes precedence)
if ~exist('inputfile','var') && ~isempty(D.gui.ini.inputfile.value)
    inputfile = D.gui.ini.inputfile.value;
end

% load the input variables
if exist('inputfile','var')
    % convert relative path to absolute
    if ~isstruct(inputfile) && ~strcmp(inputfile(1),'/') && ~strcmp(inputfile(2),':')
        inputfile = fullfile(D.gui.ini.basedir.value,inputfile);
    end
    % load the inputfile
    if isstruct(inputfile) && isfield(inputfile,'version')
        % input is the D structure
        guidata(H,inputfile);
        D = guidata(H);
        
        % update the state of the GUI
        iniset(H);drawnow
        guidata(H,D);drawnow

        updateinfo(H)
        inputfile = sprintf('var: %s',inputname(1));
    elseif exist(inputfile,'file')
        % inputfile is specified absolute
        infoload([],[],H,inputfile)
    else
        msgstr = sprintf('The input file (%s) was not found, nothing loaded',inputfile);
        msgdlgjn(msgstr,dlgposition(H));
        return
    end
else
    return
end

% get the current guidata
D = guidata(H);

% store the inputfile
D.inputfile = inputfile;

% update the status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[9] basicgdic data loaded from %s',inputfile));
guidata(H,D);drawnow

% =============================================
if ~headlessmode
    varargout{1} = D;
    return
end

% convert relative path to absolute
if ~strcmp(outputfile(1),'/') && ~strcmp(outputfile(2),':')
    outputfile = fullfile(D.gui.ini.basedir.value,outputfile);
end

D.outputfile = outputfile;
% update the application data
guidata(H,D);drawnow;


if headlessmode
    headlessstatus(sprintf('inputfile %s loaded (GUI handle is %d)',inputfile,H));
    % disable liveview (nobody is watching)
    set(S.corliveview,'Value',1)
    % start the correlation
    headlessstatus('starting correlation');
    corgo([],[],H)
    % save the result
    D = guidata(H);
    D = rmfield(D,'outputfile');
    guidata(H,D);
    headlessstatus(sprintf('writing outputfile %s',outputfile));
    infosave([],[],H,outputfile)
    headlessstatus(sprintf('outputfile %s saved',outputfile));
    % exit gui
    delete(H)
    headlessstatus('exiting basicgdic');
end

% set the output
varargout{1} = D;


% =========================================================================
% Section 0: Main GUI functions
% =========================================================================

% ==================================================
function iniwrite(filename,ini)
% this function reads all the guihandles and populates the ini structure
% from the current state of the GUI
fields = fieldnames(ini);
N = length(fields);

fid = fopen(filename,'w+t');
% print the header
fprintf(fid,'%% settings file written by basicgdic on %s \n',datestr(now));
fprintf(fid,'%% ================================================== \n');


% set all the options
for k = 1:N
    name = fields{k};
    type = ini.(fields{k}).type;
    value = ini.(fields{k}).value;
    if ~strcmpi(type,'string')
        value = num2str(value);
    end
    fprintf(fid,'%s:%s = %s\n',name,type,value);
end

fprintf(fid,'%% ================================================== \n');
fclose(fid);


% ==================================================
function ini = iniread(varargin)
filename = 'basicgdic_defaults.ini';

% create the defaults file if it doesn't exist
basicgdic_defaults;

if ~exist(filename,'file')
    ini = [];
    return
end

% read the entire option file into memory
fid = fopen(filename,'rt');
F = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
F = F{1};

% Remove empty lines
I = regexpi(F,'.*');
I = ~cellfun('isempty',I);
F = F(I);

% Remove comment lines
I = regexpi(F,'^[%#\[].*');
I = cellfun('isempty',I);
F = F(I);

% Remove inline comments
I = regexpi(F,'\s*[%#\[].*');
J = ~cellfun('isempty',I);
n = 1:length(J);
for k = n(J)
    % keep only the characters before the match
    F{k} = F{k}(1:I{k}-1);
end

% Remove extra white spaces
F = strtrim(F);

% process to a structure
for k = 1:length(F)
    % one line
    A = F{k};
    try
        % evaluate the line
        S = regexpi(A,'^\s*(.*?)\s*\:\s*(.*?)\s*=\s*(.*?)\s*$','tokens');
        name = S{1}{1};
        type = S{1}{2};
        value = S{1}{3};
    catch
        fprintf(2,'!!! iniread: error processing the line [%s]\n',A);
        return
    end
    
    % convert string to double
    if ~strcmpi(type,'string')
        value = eval(value);
    end
    
    % correct the basedir
    if strcmpi(name,'basedir')
        % change all slashes to foreward
        value = regexprep(value,'\\','/');
        if isempty(value)
            value = pwd;
        elseif ~strcmp(value(1),'/') && ~strcmp(value(2),':')
            value = fullfile(pwd,value);
        end
        
    end
    
    ini.(name).type = type;
    ini.(name).value = value;
end

% ==================================================
function ini = iniget(H)
% this function reads all the guihandles and populates the ini structure
% from the current state of the GUI
S = guihandles(H);
D = guidata(H);

% process the ini
% =============================
fields = fieldnames(D.gui.ini);
N = length(fields);
ini = D.gui.ini;

% set all the options
for k = 1:N
    if isfield(S,fields{k}) && ishandle(S.(fields{k}))
        ini.(fields{k}).value = get(S.(fields{k}),ini.(fields{k}).type);
    end
end


% ==================================================
function iniset(H)
% this function should read the contents of D.gui.ini and apply all
% properties. Additionally, it should read from D.files and D.basis to fill
% those properties (these are done first). Finally, it should correctly
% enable/disable all options and correctly set the sliders
S = guihandles(H);
D = guidata(H);

% update the Images GUI part
fileset(H)

% update the Basis GUI part
basisset(H)
 
% process the ini
% =============================
fields = fieldnames(D.gui.ini);
N = length(fields);
ini = D.gui.ini;

% set all the options
for k = 1:N
    if isfield(S,fields{k}) && ishandle(S.(fields{k}))
        set(S.(fields{k}),ini.(fields{k}).type,ini.(fields{k}).value)
    end
end

% Correct the basis choise if necessary
% =============================
basislst = get(S.basislist,'String');
Nbasis = length(basislst);
set(S.prepAbasis,'Value',min([Nbasis get(S.prepAbasis,'Value')]));
set(S.prepBbasis,'Value',min([Nbasis get(S.prepBbasis,'Value')]));
set(S.prepCbasis,'Value',min([Nbasis get(S.prepCbasis,'Value')]));
set(S.finalbasis,'Value',min([Nbasis get(S.finalbasis,'Value')]));
set(S.dicrelbasis,'Value',min([Nbasis+1 get(S.dicrelbasis,'Value')]));

% Correct the enable states
% =============================
basistype([],[],H)

cgstr{1,1} = 'prepA';
cgstr{2,1} = 'prepB';
cgstr{3,1} = 'prepC';
cgstr{4,1} = 'final';
for k = 1:3
    if ini.([cgstr{k} 'on']).value == 1
        set(S.([cgstr{k} 'on']),'BackgroundColor',D.gui.color.hl);
        set(S.([cgstr{k} 'off']),'BackgroundColor',D.gui.color.bg);
        set(S.([cgstr{k} 'on']),'Value',1);
        set(S.([cgstr{k} 'off']),'Value',0);
        enable = 'on';
    else
        set(S.([cgstr{k} 'on']),'BackgroundColor',D.gui.color.bg);
        set(S.([cgstr{k} 'off']),'BackgroundColor',D.gui.color.hl);
        set(S.([cgstr{k} 'on']),'Value',0);
        set(S.([cgstr{k} 'off']),'Value',1);
        enable = 'off';
    end
end
    
for k = 1:4
    set(S.([cgstr{k} 'blur']),'Enable',enable);
    set(S.([cgstr{k} 'level']),'Enable',enable);
    set(S.([cgstr{k} 'improcess']),'Enable',enable);
    set(S.([cgstr{k} 'basis']),'Enable',enable);
    set(S.([cgstr{k} 'convcrit']),'Enable',enable);
    set(S.([cgstr{k} 'maxit']),'Enable',enable);
    set(S.([cgstr{k} 'gradient']),'Enable',enable);
    set(S.([cgstr{k} 'tikhpar1']),'Enable',enable);
    set(S.([cgstr{k} 'tikhpar2']),'Enable',enable);
    set(S.([cgstr{k} 'tikhsteps']),'Enable',enable);
end
drawnow;

% Dimensions
if get(S.dicdimensions,'Value') == 1 %2D
    set(S.dicrelaxation,'Enable','on')
    if get(S.dicrelaxation,'Value') == 1 %none
        set(S.dicrelbasis,'Enable','off')
    else % brightness or brightness+contrast
        set(S.dicrelbasis,'Enable','on')
    end
elseif get(S.dicdimensions,'Value') == 2 %3D
    set(S.dicrelaxation,'Enable','off')
    set(S.dicrelbasis,'Enable','off')
    set(S.resstraindef,'Value',5)
end

% Relaxation
if get(S.dicrelaxation,'Value') == 1 %none
    set(S.dicrelbasis,'Enable','off')
else % brightness or brightness+contrast
    set(S.dicrelbasis,'Enable','on')
end

% Correct the sliders
% =============================
set(S.basissoftenslider,'Value',str2double(ini.basissoften.value));
set(S.basisalphaslider,'Value',str2double(ini.basisalpha.value));
set(S.diclagrangeslider,'Value',str2double(ini.diclagrange.value));
set(S.ressofteningslider,'Value',str2double(ini.ressoftening.value));
set(S.resalphaslider,'Value',str2double(ini.resalpha.value));
set(S.resarrowscaleslider,'Value',str2double(ini.resarrowscale.value));

% image sliders
sliderupdate(H)

% Correct the indicators
% =============================
if isfield(D,'files')
    set(S.secind1,'BackgroundColor',D.gui.color.on)
else
    set(S.secind1,'BackgroundColor',D.gui.color.off)
end
if isfield(D,'pateval')
    set(S.secind2,'BackgroundColor',D.gui.color.on)
else
    set(S.secind2,'BackgroundColor',D.gui.color.init)
end
if isfield(D,'roi')
    set(S.secind3,'BackgroundColor',D.gui.color.on)
else
    set(S.secind3,'BackgroundColor',D.gui.color.off)
end
if isfield(D,'iguess')
    set(S.secind4,'BackgroundColor',D.gui.color.on)
else
    set(S.secind4,'BackgroundColor',D.gui.color.init)
end
if isfield(D,'basis')
    set(S.secind5,'BackgroundColor',D.gui.color.on)
else
    set(S.secind5,'BackgroundColor',D.gui.color.off)
end
if isfield(D,'cor')
    set(S.secind6,'BackgroundColor',D.gui.color.on)
    set(S.secind7,'BackgroundColor',D.gui.color.on)
else
    set(S.secind6,'BackgroundColor',D.gui.color.init)
    set(S.secind7,'BackgroundColor',D.gui.color.off)
end
if isfield(D,'res')
    set(S.secind8,'BackgroundColor',D.gui.color.on)
else
    set(S.secind8,'BackgroundColor',D.gui.color.off)
end


% ==================================================
function fileset(H)
% this function reads D.files and populates the corresponding GUI parts

S = guihandles(H);
D = guidata(H);

if ~isfield(D,'files')
    set(S.filelist,'String',{''});
    set(S.filelist,'Max',2);
    return
end

% fix the filelist
filelist = {D.files(:).name}';
N = length(filelist);
set(S.filelist,'String',filelist);
set(S.filelist,'Max',N);

% update info
info(1).str = 'Nfiles';
info(1).val = N;
D.gui.info = appendinfo(D.gui.info,info);

guidata(H,D);drawnow


% ==================================================
function basisset(H)
% this function reads D.basis and populates the corresponding GUI parts

S = guihandles(H);
D = guidata(H);

if ~isfield(D,'basis')
    set(S.basislist,'String',{'basis 1';'basis 2';'basis 3';'basis 4'})
    set(S.basislist,'Value',1)
    return
end

% populate the basislist
Nb = length(D.basis);
done = 0;
for id = 1:Nb
    if ~isempty(D.basis(id).name)
        done = done + 1;
        basislst{id,1} = D.basis(id).name;
    else
        basislst{id,1} = sprintf('basis %d',id);
    end
end
set(S.basislist,'String',basislst)
set(S.basislist,'Value',1)

% set the enable states
basistype([],[],H)

% check if the current selection still makes sense
tags{1,1} = 'prepAbasis';
tags{2,1} = 'prepBbasis';
tags{3,1} = 'prepCbasis';
tags{4,1} = 'finalbasis';
tags{5,1} = 'dicrelbasis';
for k = 1:5
    % append to the current list
    if k == 5
        basislst = [{'same as U'} ; basislst];
    end
    % get the current string
    set(S.(tags{k}),'Value',1)
    set(S.(tags{k}),'String',basislst)
end

% Set the indicator
if done == Nb
    set(S.secind5,'BackgroundColor',D.gui.color.on)
end


% ==================================================
function sliderupdate(H)
% this function fixes the length of all sliders, such that it matches the
% number of images and basis functions.

S = guihandles(H);
D = guidata(H);

% get the number of images
if ~isfield(D,'files')
    Nfiles = 1;
else
    Nfiles = length(D.files);
end

% setting the maximum slider position
slidermax = max([Nfiles 2]);
sliderstep = [1/(slidermax-1) 10/(slidermax-1)];

% update sliders (related to images)
set(S.imageid,'String','1')
set(S.imageslider,'Value',1)
set(S.imageslider,'Min',1)
set(S.imageslider,'Max',slidermax)
set(S.imageslider,'Sliderstep',sliderstep)
set(S.patid,'String','1')
set(S.patslider,'Value',1)
set(S.patslider,'Min',1)
set(S.patslider,'Max',slidermax)
set(S.patslider,'Sliderstep',sliderstep)
set(S.roiid,'String','1')
set(S.roislider,'Value',1)
set(S.roislider,'Min',1)
set(S.roislider,'Max',slidermax)
set(S.roislider,'Sliderstep',sliderstep)
set(S.iguessid,'String','1')
set(S.iguessslider,'Value',1)
set(S.iguessslider,'Min',1)
set(S.iguessslider,'Max',slidermax)
set(S.iguessslider,'Sliderstep',sliderstep)
set(S.corid,'String','0')
set(S.corslider,'Value',0)
set(S.corslider,'Min',0)
set(S.corslider,'Max',slidermax-1)
set(S.corslider,'Sliderstep',sliderstep)
set(S.resid,'String','1')
set(S.resslider,'Value',0)
set(S.resslider,'Min',0)
set(S.resslider,'Max',slidermax-1)
set(S.resslider,'Sliderstep',sliderstep)

% get the number of basis functions
k = get(S.basislist,'Value');
if isfield(D,'basis') && isfield(D.basis,'plotphi')
    Nphi = D.basis(k).Nphi;
else
    Nphi = 1;
end

% setting the maximum slider position
slidermax = max([Nphi 2]);
if Nphi > 2
    sliderstep = [1/(slidermax-2) 10/(slidermax-2)];
else
    sliderstep = [1/(slidermax-1) 10/(slidermax-1)];
end

% update sliders (related to basis functions)
set(S.phiid,'String','1')
set(S.phislider,'Value',1)
set(S.phislider,'Max',slidermax)
set(S.phislider,'Sliderstep',sliderstep)

drawnow;



% ==================================================
function WindowScrollWheelFcn(H,evnt)
% this function is called when using the mousewheel in the gui
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end

% number of images
Nim = length(D.files);

% depending on the current panel, get the current slider value (id)
if D.gui.activepanel == 1
    id = str2double(get(S.imageid,'String'));
    id = round(id);
elseif D.gui.activepanel == 2
    id = str2double(get(S.patid,'String'));
    id = round(id);
elseif D.gui.activepanel == 3
    id = str2double(get(S.roiid,'String'));
    id = round(id);
elseif D.gui.activepanel == 4
    id = str2double(get(S.iguessid,'String'));
    id = round(id);
elseif D.gui.activepanel == 5
    id = str2double(get(S.phiid,'String'));
    id = round(id);
    k = get(S.basislist,'Value');
    if isfield(D,'basis') && isfield(D.basis(k),'Nphi')
        Nim = D.basis(k).Nphi;
        Nim = max([Nim 2]);
    else
        return
    end
elseif D.gui.activepanel == 7
    id = str2double(get(S.corid,'String'));
    id = round(id);
elseif D.gui.activepanel == 8
    id = str2double(get(S.resid,'String'));
    id = round(id);
else
    return
end

% mouse wheel up or down
if evnt.VerticalScrollCount < 0
    id = id - 1 ;
elseif evnt.VerticalScrollCount > 0
    id = id + 1 ;
end

% limit the id to the slider limits
if any(D.gui.activepanel == [1 2 3 4 5 6 9])
    id = max([id 1]);
    id = min([id Nim]);
elseif any(D.gui.activepanel == [7,8])
    id = max([id 0]);
    id = min([id Nim-1]);
end

% depending on the current panel, change the correct slider
if D.gui.activepanel == 1
    set(S.imageslider,'Value',id);
    set(S.imageid,'String',num2str(id));
    plotimage(H)
elseif D.gui.activepanel == 2
    set(S.patslider,'Value',id);
    set(S.patid,'String',num2str(id));
    if D.gui.activepatview == 1
        plotpattern(H)
    elseif D.gui.activepatview == 2
        plotacf(H)
    end
elseif D.gui.activepanel == 3
    set(S.roislider,'Value',id);
    set(S.roiid,'String',num2str(id));
    plotroimask(H)
elseif D.gui.activepanel == 4
    set(S.iguessslider,'Value',id);
    set(S.iguessid,'String',num2str(id));
    plotiguess(H)
elseif D.gui.activepanel == 5
    set(S.phislider,'Value',id);
    set(S.phiid,'String',num2str(id));
    plotbasis([],[],H);
elseif D.gui.activepanel == 6
elseif D.gui.activepanel == 7
    set(S.corslider,'Value',id);
    set(S.corid,'String',num2str(id));
    plotcor(H)
elseif D.gui.activepanel == 8
    set(S.resslider,'Value',id);
    set(S.resid,'String',num2str(id));
    plotres([],[],H)
end
drawnow

% ==================================================
function unlockGUI(H,evnt)
% function called when pressing the unlock toolbar button
disp('unlocking')
D = guidata(H);
set(D.gui.waithandles,'Enable','on');drawnow


% ==================================================
function keyPressFcn(H,evnt)
% call function for keypresses in the gui
D = guidata(H);

% press escale to "exit" a blocked state
if strcmp(evnt.Key,'escape')
    % enable gui controls
    disp('unlocking')
    set(D.gui.waithandles,'Enable','on');drawnow
end

% Ctrl-N shortcuts to the sections
if strcmp(evnt.Modifier,'control')
    if strcmp(evnt.Key,'1')
        secbutton([],[],H,1)
    elseif strcmp(evnt.Key,'2')
        secbutton([],[],H,2)
    elseif strcmp(evnt.Key,'3')
        secbutton([],[],H,3)
    elseif strcmp(evnt.Key,'4')
        secbutton([],[],H,4)
    elseif strcmp(evnt.Key,'5')
        secbutton([],[],H,5)
    elseif strcmp(evnt.Key,'6')
        secbutton([],[],H,6)
    elseif strcmp(evnt.Key,'7')
        secbutton([],[],H,7)
    elseif strcmp(evnt.Key,'8')
        secbutton([],[],H,8)
    elseif strcmp(evnt.Key,'9')
        secbutton([],[],H,9)
    end
end

% arrows
if strcmp(evnt.Key,'leftarrow') || strcmp(evnt.Key,'downarrow')
    E.VerticalScrollCount = -1;
    WindowScrollWheelFcn(H,E)
elseif strcmp(evnt.Key,'rightarrow') || strcmp(evnt.Key,'uparrow')
    E.VerticalScrollCount = 1;
    WindowScrollWheelFcn(H,E)
end

% ==================================================
function [] = secbutton(varargin)
% Button Callback: Switch section pannels
H = varargin{3};
S = guihandles(H);
D = guidata(H);
% the pressed button
but = varargin{4};
D.gui.activepanel = but;

% switch all buttons to gray, and one to blue
for k = 1:9
    if k == but
        set(S.(['secpan' num2str(k)]),'Visible','on');
        set(S.(['figpan' num2str(k)]),'Visible','on');
        set(S.(['secbut' num2str(k)]),'BackgroundColor',D.gui.color.hl);
    else
        set(S.(['secpan' num2str(k)]),'Visible','off');
        set(S.(['figpan' num2str(k)]),'Visible','off');
        set(S.(['secbut' num2str(k)]),'BackgroundColor',D.gui.color.bg);
    end
end

% change colormap between gray and colored
if any(but == [5 8])
    colormap(D.gui.color.cmapoverlay)
else
    colormap(D.gui.color.cmap)
end

% set indicator on for DIC options
if but == 6
    set(S.secind6,'BackgroundColor',D.gui.color.on)
elseif but == 9
    set(S.secind9,'BackgroundColor',D.gui.color.on)
end

% update the application data
guidata(H,D);

% update the corresponding figures
if D.gui.activepanel == 1
    cla(S.axes1);
    plotimage(H);
elseif D.gui.activepanel == 2
    cla(S.axes2);
    if D.gui.activepatview == 1
        plotpattern(H)
    elseif D.gui.activepatview == 2
        plotacf(H)
    end
elseif D.gui.activepanel == 3
    cla(S.axes3);
    plotroimask(H);
elseif D.gui.activepanel == 4
    cla(S.axes41);
    cla(S.axes42);
    cla(S.axes43);
    cla(S.axes44);
    plotiguess(H);
elseif D.gui.activepanel == 5
    cla(S.axes4);
    plotbasis([],[],H);
elseif D.gui.activepanel == 6
    % update the basis list
    basislst = get(S.basislist,'String');
    % check if the current selection still makes sense
    tags{1,1} = 'prepAbasis';
    tags{2,1} = 'prepBbasis';
    tags{3,1} = 'prepCbasis';
    tags{4,1} = 'finalbasis';
    tags{5,1} = 'dicrelbasis';
    for k = 1:5
        % get the current string
        strs = get(S.(tags{k}),'String');
        val  = get(S.(tags{k}),'Value');
        val = min([val length(strs)]);
        val = max([val 1]);
        str = strs{val};
        % append to the current list
        if k == 5
            basislst = [{'same as U'} ; basislst];
        end
        id = find(strcmp(str,basislst),1);
        if isempty(id)
            val = min([val Nb]);
            val = max([val 1]);
        else
            val = id;
        end
        set(S.(tags{k}),'Value',val)
        set(S.(tags{k}),'String',basislst)
    end

elseif D.gui.activepanel == 7
    cla(S.axes7);
    plotcor(H);
elseif D.gui.activepanel == 8
    cla(S.axes8);
    plotres([],[],H);
elseif D.gui.activepanel == 9
    updateinfo(H)
end

% force redrawing of gui elements
drawnow


% ==================================================
function [] = infobutton(varargin)
% Button Callback: Switch info panels (Info,Status,Help)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
but = varargin{4};
D.gui.activeinfopanel = but;

% switch all buttons to gray, and one to blue
for k = 1:3
    if k == but
        set(S.(['infopan' num2str(k)]),'Visible','on');
        set(S.(['infobut' num2str(k)]),'BackgroundColor',D.gui.color.hl);
    else
        set(S.(['infopan' num2str(k)]),'Visible','off');
        set(S.(['infobut' num2str(k)]),'BackgroundColor',D.gui.color.bg);
    end
end

% update the application data
guidata(H,D);




% =========================================================================
% Section 1: Images
% =========================================================================

% ==================================================
function [] = fileadd(varargin)
% Load additional images

H = varargin{3};
S = guihandles(H);
D = guidata(H);

% old filelist
filelist = get(S.filelist,'String');

% old file datastructure
if (isfield(D,'files')) && (~isempty(D.files))
    files = D.files;
    Nfiles = length(files);
    val = get(S.filelist,'Value');
    if ~isempty(val)
        basedir = files(val(end)).path;
    else
        basedir = files(end).path;
    end
else
    Nfiles = 0;
    % get the previous basedir
    if ~isempty(D.gui.ini.basedir.value);
        basedir = D.gui.ini.basedir.value;
    else
        basedir = pwd;
    end
end

% get the filesize of the reference
if isfield(D,'files')
    nref = D.files(1).size(1);
    mref = D.files(1).size(2);
end

%uigetfile
[FileName,PathName] = uigetfile(...
    {'*.bmp;*.jpg;*.tif;*.png;*.gif;*.plu;*.BMP;*.JPG;*.TIF;*.PNG;*.GIF;*.PLU',...
    'All Image Files';...
    '*.*','All Files' },...
    'add files',...
    'MultiSelect','on',...
    basedir);

% if cancel
if isequal(FileName,0)
    return
end

% for output of uigetfile to be a string, even if only one file
if ~iscell(FileName)
    fname{1} = FileName;
else
    fname = FileName;
end

% number of new files
Nselected = length(fname);

colorimage = false;

for isel = 1:Nselected
    id = Nfiles+isel;
    bcwaitbar(H,isel/Nselected,sprintf('loading files (%d/%d)',isel,Nselected));
    
    % add files to list and datastruct
    ext = regexp(fname{isel},'\..*$','match');
    ext = ext{1}(2:end);
    files(id).name = fname{isel};
    files(id).path = PathName;
    files(id).ext = ext;
    filelist{id} = fname{isel};
    
    % read the images
    filename = fullfile(files(id).path,files(id).name);
    
    % get the extention
    ext = regexp(filename,'\.','split');
    ext = ext{end};
    
    if strcmp(ext,'plu')
        plu = pluread(filename);
        A = plu.Z;% + plu.z0;
        [n, m] = size(A);
        cls = 'double';
        
        % check for NaN values
        if any(isnan(A(:)));
            [X, Y] = meshgrid(1:m,1:n);
            I = find(~isnan(A));
            A = griddata(X(I),Y(I),A(I),X,Y);
        end
        
        set(S.respixelsize,'String',num2str(plu.pixelsize(1)));
        set(S.resunit,'Value',2);
    else
        % read image
        A = imread(filename);
        % store current class
        cls = class(A);
        % convert to double
        A = double(A);
        
        % convert RGB to Grayscale
        [n, m, k] = size(A);
        if k == 3
            colorimage = true;
            A = 0.2989*A(:,:,1) + 0.5870*A(:,:,2) + 0.1140*A(:,:,3);
        end
    end
    
    % get the image size of the reference image
    if id == 1
        nref = n;
        mref = m;
    else
        if (n ~= nref) || (m ~= mref)
            msgstr{1,1} = sprintf('Incompatible image size of %s (%d,%d), expected (%d,%d)',fname{isel},n,m,nref,mref);
            msgstr{2,1} = 'Add file procedure stopped';
            msgdlgjn(msgstr,dlgposition(H));
            delete(hwait)
            return
        end
    end
    
    % store images
    files(id).size = size(A);
    files(id).range = [min(A(:)) max(A(:))];
    files(id).class = cls;
    files(id).image = A;
    
    
    % add zero initial guess
    if isfield(D,'iguess')
        D.iguess.ux(:,id) = 0;
        D.iguess.uy(:,id) = 0;
    end
    
    % add initial guess (processed) image
    if isfield(files,'imageproc')
        % get settings
        Rblur = str2double(get(S.iguessblur,'String'));
        gradtype = get(S.iguessimprocesstype,'String');
        gradtype = gradtype(get(S.iguessimprocesstype,'Value'));
        
        % blur
        % ===========
        if Rblur > 0.3
            Im = image_blur(A,Rblur);
        end
        
        % gradient transform
        % ===========
        if ~strcmp(gradtype,'none')
            [dImdx, dImdy] = gradient(Im);
        end
        
        if strcmp(gradtype,'grad|xy|');
            Im = sqrt(dImdx.^2+dImdy.^2);
        elseif strcmp(gradtype,'grad|x|');
            Im = abs(dImdx);
        elseif strcmp(gradtype,'grad|y|');
            Im = abs(dImdy);
        elseif strcmp(gradtype,'grad(xy)');
            Im = 0.5*(dImdx+dImdy);
        elseif strcmp(gradtype,'grad(x)');
            Im = dImdx;
        elseif strcmp(gradtype,'grad(y)');
            Im = dImdy;
        end
        
        files(id).imageproc = Im;
    end

    % add zero initial guess
    if isfield(D,'cor')
        D.cor(id-1).inc = id-1;
    end
    
    % add zero initial guess
    if isfield(D,'res')
        D.res(id-1).inc = id-1;
    end
    
end
bcwaitbar(H);


if colorimage
    msgstr{1,1} = 'Color images detected and converted to grayscale.';
    msgstr{2,1} = 'The blur radii of the coarse grain steps are set to 2 pixels';
    msgstr{3,1} = 'to reduce demosaicing artifacts';
    msgdlgjn(msgstr,dlgposition(H));
    
    tags{1,1} = 'prepAblur';
    tags{2,1} = 'prepBblur';
    tags{3,1} = 'prepCblur';
    tags{4,1} = 'finalblur';
    for k = 1:4
        set(S.(tags{k}),'String','2')
    end
end

% update list in gui
Nfiles = length(files);
set(S.filelist,'String',filelist);
set(S.filelist,'Max',Nfiles);
D.files = files;

% update info
info(1).str = 'Nfiles';
info(1).val = Nfiles;
info(2).str = 'n';
info(2).val = n;
info(3).str = 'm';
info(3).val = m;
info(4).str = 'Npx';
info(4).val = n*m;
info(5).str = 'class';
info(5).val = cls;
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[1] %d files added',Nselected));

% update the application data
guidata(H,D);

% update sliders
sliderupdate(H);

% plot the figures
plotimage(H);

% set section indicator on
set(S.secind1,'BackgroundColor',D.gui.color.on)



% ==================================================
function [] = filedel(varargin)
% remove files from the list
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if isfield(D,'files')
    files = D.files;
    Nfiles = length(files);
else
    return
end

% get current filelist
filelist = get(S.filelist,'String');
% get selection
selected = get(S.filelist,'Value');
Nselected = length(selected);

if Nselected == 1
    filestr = 'file';
else
    filestr = 'files';
end

% ask for confirmation
if any(selected==1)
    qstr{1,1} = 'Warning: Removing the reference image (first file) will reset all results';
    qstr{2,1} = sprintf('Delete %d %s from the list?',Nselected,filestr);
else
    qstr{1,1} = sprintf('Delete %d %s from the list?',Nselected,filestr);
end
button = questdlgjn(qstr,'deleting files','Yes','No','Yes',dlgposition(H));
if isempty(button) || strcmpi(button,'No')
    return
end

if any(selected==1)
    set(S.secind1,'BackGroundColor',D.gui.color.off)
    set(S.secind2,'BackGroundColor',D.gui.color.init)
    set(S.secind3,'BackGroundColor',D.gui.color.off)
    set(S.secind4,'BackGroundColor',D.gui.color.init)
    set(S.secind5,'BackGroundColor',D.gui.color.off)
    set(S.secind6,'BackGroundColor',D.gui.color.init)
    set(S.secind7,'BackGroundColor',D.gui.color.off)
    set(S.secind8,'BackGroundColor',D.gui.color.init)
    set(S.secind9,'BackGroundColor',D.gui.color.init)
end


if Nselected == Nfiles
    % delete all files
    D = rmfield(D,'files');
    Nfiles = 2;
    filelist = {};
    
    if isfield(D,'iguess')
        D = rmfield(D,'iguess');
    end
    if isfield(D,'cor')
        D = rmfield(D,'cor');
    end
    if isfield(D,'res')
        D = rmfield(D,'res');
    end
    if isfield(D,'roi')
        D = rmfield(D,'roi');
    end
    if isfield(D,'pateval')
        D = rmfield(D,'pateval');
    end
else
    % invert index
    I = 1:Nfiles;
    I = setdiff(I,selected);
    
    Inc = 1:(Nfiles-1);
    Inc = setdiff(Inc,selected-1);
    
    % delete files from lists
    filelist = filelist(I);
    files = files(I);
    Nfiles = length(files);
    D.files = files;
    % Delete corresponding data
    if isfield(D,'iguess')
        if any(selected == 1) % reference is deleted, start over
            D = rmfield(D,'iguess');
        else
            D.iguess.ux = D.iguess.ux(:,I);
            D.iguess.uy = D.iguess.uy(:,I);
        end
    end
    if isfield(D,'cor')
        if any(selected == 1) % reference is deleted, start over
            D = rmfield(D,'cor');
        else
            D.cor = D.cor(Inc);
        end
    end
    if isfield(D,'res')
        if any(selected == 1) % reference is deleted, start over
            D = rmfield(D,'res');
        else
            D.res = D.res(Inc);
        end
    end
end

% update gui
set(S.filelist,'Value',1);
set(S.filelist,'Max',Nfiles);
set(S.filelist,'String',filelist);

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[1] %d files deleted',Nselected));

% update the application data
guidata(H,D);

% update sliders
sliderupdate(H);

% plot the figures
plotimage(H);

% ==================================================
function [] = fileup(varargin)
% move files up the list
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if isfield(D,'files')
    files = D.files;
    Nfiles = length(files);
else
    return
end

% get current filelist
filelist = get(S.filelist,'String');
% get selection
sel = get(S.filelist,'Value');
Nsel = length(sel);

% files above selected
ab = sel-1;
ab(ab<=1) = 1;

% index
I = 1:Nfiles;
Inc = I - 1;

% swap selected index with item above
for k = 1:length(sel)
    I([sel(k),ab(k)])   = I([ab(k),sel(k)]);
    Inc([sel(k),ab(k)]) = Inc([ab(k),sel(k)]);
end
Inc = Inc(2:end);

if Nsel == 1
    filestr = 'file';
else
    filestr = 'files';
end

% ask for confirmation
if I(1) ~= 1 % changing the reference
    qstr{1,1} = 'Warning: Changing the reference image (first file) will reset all results';
    qstr{2,1} = sprintf('Move %d %s in list?',Nsel,filestr);
    button = questdlgjn(qstr,'move files','Yes','No','Yes',dlgposition(H));
    if isempty(button) || strcmpi(button,'No')
        return
    end
end

% update to new list order
filelist = filelist(I);
files = files(I);
Nfiles = length(files);

% update gui
set(S.filelist,'Value',ab);
set(S.filelist,'String',filelist);
D.files = files;

% Move corresponding data
if isfield(D,'iguess')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'iguess');
    else
        D.iguess.ux = D.iguess.ux(:,I);
        D.iguess.uy = D.iguess.uy(:,I);
    end
end
if isfield(D,'cor')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'cor');
    else
        D.cor = D.cor(Inc);
    end
end
if isfield(D,'res')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'res');
    else
        D.res = D.res(Inc);
    end
end

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[1] %d files moved up in the list',Nsel));

% update the application data
guidata(H,D);

% plot the figures
plotimage(H);


% ==================================================
function [] = filedown(varargin)
% move files down the list

H = varargin{3};
S = guihandles(H);
D = guidata(H);

if isfield(D,'files')
    files = D.files;
    Nfiles = length(files);
else
    return;
end

% get current filelist
filelist = get(S.filelist,'String');
% get selection
sel = get(S.filelist,'Value');
Nsel = length(sel);

% file below selected
bel = sel+1;
bel(bel>=Nfiles) = Nfiles;

% index
I = 1:Nfiles;
Inc = I-1;

% switch index with item below
for k = length(sel):-1:1
    I([sel(k),bel(k)])   = I([bel(k),sel(k)]);
    Inc([sel(k),bel(k)]) = Inc([bel(k),sel(k)]);
end
Inc = Inc(2:end);

if Nsel == 1
    filestr = 'file';
else
    filestr = 'files';
end

% ask for confirmation
if I(1) ~= 1 % changing the reference
    qstr{1,1} = 'Warning: Changing the reference image (first file) will reset all results';
    qstr{2,1} = sprintf('Move %d %s in list?',Nsel,filestr);
    button = questdlgjn(qstr,'move files','Yes','No','Yes',dlgposition(H));
    if isempty(button) || strcmpi(button,'No')
        return
    end
end

% update file order
filelist = filelist(I);
files = files(I);

% update gui
set(S.filelist,'Value',bel);
set(S.filelist,'String',filelist);
D.files = files;

% Move corresponding data
if isfield(D,'iguess')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'iguess');
    else
        D.iguess.ux = D.iguess.ux(:,I);
        D.iguess.uy = D.iguess.uy(:,I);
    end
end
if isfield(D,'cor')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'cor');
    else
        D.cor = D.cor(Inc);
    end
end
if isfield(D,'res')
    if I(1) ~= 1 % changing the reference
        D = rmfield(D,'res');
    else
        D.res = D.res(Inc);
    end
end

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[1] %d files moved down in the list',Nsel));

% update the application data
guidata(H,D);

% plot the figures
plotimage(H);

% ==================================================
function [] = imageid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.imageid,'String'));
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.imageslider,'Value',id);
plotimage(H)

% ==================================================
function [] = fileselect(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.filelist,'Value');
id = id(1);
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.imageid,'String',num2str(id));
set(S.imageslider,'Value',id);
plotimage(H)


% ==================================================
function [] = imageslider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.imageslider,'Value');
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.imageid,'String',num2str(id));
plotimage(H)


% =========================================================================
% Section 2: Pattern Evaluation
% =========================================================================

% ==================================================
function [] = patid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.patid,'String'));
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.patslider,'Value',id);
plotpattern(H)

% ==================================================
function [] = patslider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.patslider,'Value');
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.patid,'String',num2str(id));
plotpattern(H)

% ==================================================
function [] = patshowpat(varargin)
% Button Callback: Show the pattern
H = varargin{3};
S = guihandles(H);
D = guidata(H);
set(S.patshowpat,'BackgroundColor',D.gui.color.hl);
set(S.patshowACF,'BackgroundColor',D.gui.color.bg);
plotpattern(H);
D.gui.activepatview = 1;
guidata(H,D);

% ==================================================
function [] = patshowACF(varargin)
% Button Callback: Show the autocorrelation function
H = varargin{3};
S = guihandles(H);
D = guidata(H);
set(S.patshowpat,'BackgroundColor',D.gui.color.bg);
set(S.patshowACF,'BackgroundColor',D.gui.color.hl);
plotacf(H);
D.gui.activepatview = 2;
guidata(H,D);

% =========================================================================
% Section 3: ROI and Mask
% =========================================================================

% ==================================================
function [] = roiid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.roiid,'String'));
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.roislider,'Value',id);
plotroimask(H)

% ==================================================
function [] = roislider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.roislider,'Value');
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.roiid,'String',num2str(id));
plotroimask(H)

% ==================================================
function [] = masksave(varargin)
% Button Callback: save the mask
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uiputfile(fullfile(basedir,'mask.mat'), 'save the mask');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;


% save to file
bgdic.mask     = D.mask;
bgdic.roi      = D.roi;
bgdic.version  = D.version;
bgdic.savetype = 'mask';
save('-v7.3',fullfile(pathname,filename),'bgdic')

% update status
D.gui.stat = appendstatus(D.gui.stat,'[3] Mask saved to file');
guidata(H,D);

set(D.gui.waithandles,'Enable','on');drawnow



% ==================================================
function [] = maskload(varargin)
% Button Callback: load the mask
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uigetfile(fullfile(basedir,'mask.mat'), 'load the mask');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

% load from file
load(fullfile(pathname,filename),'bgdic')
if ~exist('bgdic','var') || ~isfield(bgdic,'savetype')
    msgstr = 'incompatible save file';
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end
if ~strcmp(bgdic.savetype,'mask')
    msgstr = sprintf('incorrect save type (%s), try loading it in another section',bgdic.savetype);
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

D.mask = bgdic.mask;
D.roi = bgdic.roi;

if isfield(D,'savetype')
    D = rmfield(D,'savetype');
end

% update status
D.gui.stat = appendstatus(D.gui.stat,'[3] Mask loaded from file');

% update the application data
guidata(H,D);

set(D.gui.waithandles,'Enable','on');drawnow

plotroimask(H);


% ==================================================
function [] = roiset(varargin)
% Button Callback: Define a region of interest
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
if isfield(D,'cor')
    qstr{1,1} = 'Warning: Changing the ROI will reset all results';
    qstr{2,1} = 'Continue?';
    button = questdlgjn(qstr,'Change the ROI','Yes','No','Yes',dlgposition(H));
    if isempty(button) || strcmpi(button,'No')
        return
    end
    D = rmfield(D,'cor');
end
if isfield(D,'res')
    D = rmfield(D,'res');
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

n = D.files(1).size(1);
m = D.files(1).size(2);
% Select an area
set(H,'CurrentAxes',S.axes3);
h = title('position the rectangle, confirm with a doubleclick');
set(h,'color',D.gui.color.axeshl);

% reuse roi
if isfield(D,'roi')
    D.roi = max([1   2 1   2;D.roi]);
    D.roi = min([m-1 m n-1 n;D.roi]);
    rect = [D.roi(1),D.roi(3) ; D.roi(2),D.roi(4)];
else
    rect = [0.1*m,0.1*n ; 0.9*m,0.9*n];
end

% load interactive rectangle tool
position = selectarea(rect);
% reset the title
set(h,'color','k');
title('');

roi(1) = min(position(:,1));
roi(2) = max(position(:,1));
roi(3) = min(position(:,2));
roi(4) = max(position(:,2));
D.roi = roi;

% update info
info(1).str = 'roi1 (x1)';
info(1).val = roi(1);
info(2).str = 'roi2 (x2)';
info(2).val = roi(2);
info(3).str = 'roi3 (y1)';
info(3).val = roi(3);
info(4).str = 'roi4 (y2)';
info(4).val = roi(4);
D.gui.info = appendinfo(D.gui.info,info);

% set section indicator on
set(S.secind3,'BackgroundColor',D.gui.color.on)


% disable gui controls
set(D.gui.waithandles,'Enable','on');drawnow

% update status
D.gui.stat = appendstatus(D.gui.stat,'[3] Mask saved to file');

% update the application data
guidata(H,D);

plotroimask(H);

% ==================================================
function [] = maskset(varargin)
% Button Callback: Add pixels to the Mask
H = varargin{3};
modtype = varargin{4};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% load the current mask
mask = D.mask;

% image size
n = D.files(1).size(1);
m = D.files(1).size(2);

% coordinate vectors
x = 1:m;
y = 1:n;
[X, Y] = meshgrid(x,y);

% get the current tool type
tool = get(S.masktool,'String');
tool = tool{get(S.masktool,'Value')};

if strcmpi(tool,'rect')
    % update figure title
    set(H,'CurrentAxes',S.axes3)
    title('position the rectangle, confirm with a doubleclick');
    
    % load interactive rectangle tool
    B = [0.4*m,0.4*n ; 0.2*m,0.2*n];
    A = selectarea(B);
    % reset the title
    title('');
    
    xlim(1) = min(A(:,1));
    xlim(2) = max(A(:,1));
    ylim(1) = min(A(:,2));
    ylim(2) = max(A(:,2));
    
    % get the indices inside the rect
    I = find(X > xlim(1) & X < xlim(2) & Y > ylim(1) & Y < ylim(2));
    
elseif strcmpi(tool,'ellipse')
    % update figure title
    set(H,'CurrentAxes',S.axes3)
    title('position the ellipse, confirm with a doubleclick');
    
    % load interactive ellipse tool
    B = [0.3*m,0.3*n ; 0.3*m,0.1*n ; 0.25*m,0.3*n];
    [A, P] = selectarea(B,'ellipse');
    % reset the title
    title('');
    
    I = inpolygon(X,Y,P(:,1),P(:,2));
    I = find(I==1);
elseif strcmpi(tool,'polygon')
    % update figure title
    set(H,'CurrentAxes',S.axes3)
    title({'draw the polygon, confirm with a doubleclick';'(right-click to add contol points)'});
    
    % load interactive polygon tool
    B = [0.2*m,0.2*n ; 0.4*m,0.2*n ; 0.4*m,0.4*n ; 0.2*m,0.4*n ; 0.2*m,0.2*n];
    A = selectarea(B,'Polyline');
    % reset the title
    title('');
    
    I = inpolygon(X,Y,A(:,1),A(:,2));
    I = find(I==1);
elseif strcmpi(tool,'bspline')
    % update figure title
    set(H,'CurrentAxes',S.axes3)
    title({'draw the bspline, confirm with a doubleclick';'(right-click to add contol points)'});
    
    % load interactive polygon tool
    B = [0.2*m,0.2*n ; 0.4*m,0.2*n ; 0.4*m,0.4*n ; 0.2*m,0.4*n ; 0.2*m,0.2*n];
    [A, P] = selectarea(B,'BSpline',3);
    % reset the title
    title('');

    I = inpolygon(X,Y,P(:,1),P(:,2));
    I = find(I==1);
elseif strncmpi(tool,'spot',4)
    % update figure title
    set(H,'CurrentAxes',S.axes3)
    title('select a spot');
    
    % ask the user to select a spot
    B = ginputjn(1);
    % reset the title
    title('');
    
    if strcmpi(tool,'spot 5px')
        R = 5;
    elseif  strcmpi(tool,'spot 10px')
        R = 10;
    elseif  strcmpi(tool,'spot 50px')
        R = 50;
    end

    I = find( (X - B(1,1)).^2 + (Y - B(1,2)).^2 < R^2);
end

if strcmp(modtype,'add')
    mask = union(mask(:),I(:));
elseif strcmp(modtype,'del')
    mask = setdiff(mask(:),I(:));
elseif strcmp(modtype,'intersect')
    mask = intersect(mask(:),I(:));
end
mask = mask(:);

% maskimage
Im = maskimage(mask,n,m);

D.mask = mask;

% update info
info(1).str = 'mask Npx';
info(1).val = length(mask);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,['[3] Mask modified (' modtype ')']);

% update the application data
guidata(H,D);

set(D.gui.waithandles,'Enable','on');drawnow

plotroimask(H);

% ==================================================
function [] = maskinvert(varargin)
% Button Callback: Invert the mask
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% image size
n = D.files(1).size(1);
m = D.files(1).size(2);

% get the old mask
mask = D.mask;

% maskimage
Im = maskimage(mask,n,m);

% invert
Im = ~Im;

% back to mask index
D.mask = maskindex(Im);

% update info
info(1).str = 'mask Npx';
info(1).val = length(D.mask);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[3] Mask modified (invert)');

% update the application data
guidata(H,D);

plotroimask(H);

% ==================================================
function [] = maskclear(varargin)
% Button Callback: clear the mask
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end

% back to mask index
D.mask = [];

% update info
info(1).str = 'mask Npx';
info(1).val = length(D.mask);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[3] Mask cleared');

% update the application data
guidata(H,D);

plotroimask(H);

% =========================================================================
% Section 4: Initial Guess
% =========================================================================

% ==================================================
function [] = iguesssave(varargin)
% Button Callback: save the initial guess
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess,';'create an initial guess in section 4'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uiputfile(fullfile(basedir,'iguess.mat'), 'save the initial guess');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

% save to file
bgdic.iguess   = D.iguess;
bgdic.version  = D.version;
bgdic.savetype = 'iguess';
save('-v7.3',fullfile(pathname,filename),'bgdic')

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Initial guess saved to file');

% update the application data
guidata(H,D);


set(D.gui.waithandles,'Enable','on');drawnow


% ==================================================
function [] = iguessload(varargin)
% Button Callback: load the initial guess
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uigetfile(fullfile(basedir,'iguess.mat'), 'save the initial guess');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

% load from file
load(fullfile(pathname,filename),'bgdic')
if ~exist('bgdic','var') || ~isfield(bgdic,'savetype')
    msgstr = 'incompatible save file';
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end
if ~strcmp(bgdic.savetype,'iguess')
    msgstr = sprintf('incorrect save type (%s), try loading it in another section',bgdic.savetype);
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

D.iguess = bgdic.iguess;

if isfield(D,'savetype')
    D = rmfield(D,'savetype');
end


% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Initial guess loaded from file');

% update the application data
guidata(H,D);

set(D.gui.waithandles,'Enable','on');drawnow

plotiguess(H);

% ==================================================
function [] = iguessid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.iguessid,'String'));
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.iguessslider,'Value',id);
plotiguess(H)

% ==================================================
function [] = iguessslider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.iguessslider,'Value');
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nim]);
% update the gui
set(S.iguessid,'String',num2str(id));
plotiguess(H)

% ==================================================
function [] = iguessimprocess(varargin)
% Button Callback: process the images
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Image processing started');

N = length(D.files);
n = D.files(1).size(1);
m = D.files(1).size(2);
x = 1:m;
y = 1:n;

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get settings
Rblur = str2double(get(S.iguessblur,'String'));
gradtype = get(S.iguessimprocesstype,'String');
gradtype = gradtype(get(S.iguessimprocesstype,'Value'));


for k = 1:N
    bcwaitbar(H,k/N,sprintf('processing images (%d/%d)',k,N));
    Im = D.files(k).image;
    % blur
    % ===========
    if Rblur > 0.3
        Im = image_blur(Im,Rblur);
    end
    
    % gradient transform
    % ===========
    if ~strcmp(gradtype,'none')
        [dImdx, dImdy] = gradient(Im);
    end
    
    if strcmp(gradtype,'grad|xy|');
        Im = sqrt(dImdx.^2+dImdy.^2);
    elseif strcmp(gradtype,'grad|x|');
        Im = abs(dImdx);
    elseif strcmp(gradtype,'grad|y|');
        Im = abs(dImdy);
    elseif strcmp(gradtype,'grad(xy)');
        Im = 0.5*(dImdx+dImdy);
    elseif strcmp(gradtype,'grad(x)');
        Im = dImdx;
    elseif strcmp(gradtype,'grad(y)');
        Im = dImdy;
    end
    
    % store
    D.files(k).imageproc = Im;
end


% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Image processing done');

% update the application data
guidata(H,D);

% enable gui controls
set(D.gui.waithandles,'Enable','on');
bcwaitbar(H);

plotiguess(H);


% ==================================================
function [] = iguessaddauto(varargin)
% Button Callback: apply local DIC to get initial guess
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

if isfield(D.files,'imageproc')
    Img = {D.files.imageproc};
elseif isfield(D.files,'image')
    Img = {D.files.image};
else
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Auto initial guess routine started');

N = length(D.files);
n = D.files(1).size(1);
m = D.files(1).size(2);
x = 1:m;
y = 1:n;

% get previous guesses
if isfield(D,'iguess')
    iguess = D.iguess;
    Niguess = length(iguess.x);
else
    iguess.x = [];
    iguess.y = [];
    iguess.ux = [];
    iguess.uy = [];
    iguess.cc = [];
    Niguess = 0;
end

% % disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% Ask for iguess options
dlg_title = 'Automatic initial guess options';
prompt = {'Subset Size',...
          'Search Size',...
          'Subset Rows',...
          'Subset Cols',...
          'Updated'};
def{1} = D.gui.ini.iguesssubsetsize.value;
def{2} = D.gui.ini.iguesssearchsize.value;
def{3} = D.gui.ini.iguesssubsetrows.value;
def{4} = D.gui.ini.iguesssubsetcols.value;
def{5} = D.gui.ini.iguessupdated.value;
num_lines = 1;
answer = inputdlgjn(prompt,dlg_title,num_lines,def,dlgposition(H));
if isempty(answer)
    set(D.gui.waithandles,'Enable','on');
    return
end
Lzoi = eval(answer{1});
Lsearch = eval(answer{2});
Nrows = eval(answer{3});
Ncols = eval(answer{4});
updated = eval(answer{5});

% store the defaults
D.gui.ini.iguesssubsetsize.value = answer{1};
D.gui.ini.iguesssearchsize.value = answer{2};
D.gui.ini.iguesssubsetrows.value = answer{3};
D.gui.ini.iguesssubsetcols.value = answer{4};
D.gui.ini.iguessupdated.value = answer{5};

% re-use roi
if isfield(iguess,'box')
    box = iguess.box;
    rect = [box(1),box(3) ; box(2),box(4)];
elseif isfield(D,'roi')
    rect = [D.roi(1),D.roi(3) ; D.roi(2),D.roi(4)];
else
    rect = [0.1*m 0.1*n 0.8*m 0.8*n];
end

% Select an area
set(H,'CurrentAxes',S.axes41);
h = title('position the rectangle, confirm with a doubleclick');
set(h,'color',D.gui.color.axeshl);

% load interactive rectangle tool
position = selectarea(rect);

% reset the title
set(h,'color','k');
title('');

% store the rectangle for later use
box(1) = min(position(:,1));
box(2) = max(position(:,1));
box(3) = min(position(:,2));
box(4) = max(position(:,2));
iguess.box = box;

% ZOI must be at least 10 pixels
if Lzoi < 10
    Lzoi = 10;
end
% Search ZOI must be bigger than ZOI
if Lsearch < round(1.1*Lzoi)
    Lsearch = round(1.1*Lzoi);
end

Nrows = max([Nrows 1]);
Ncols = max([Ncols 1]);
Nzoi = Ncols*Nrows;

% create subset center locations
boxedge = min([D.gui.ini.iguessboxedge.value 0.1*Lzoi]);
if Nzoi == 1
    xc = mean(box([1,2]));
    yc = mean(box([3,4]));
elseif Nrows == 1
    xc = linspace(box(1)+boxedge,box(2)-boxedge,Ncols);
    yc = mean(box([3,4]));
elseif Ncols == 1
    xc = mean(box([1,2]));
    yc = linspace(box(3)+boxedge,box(4)-boxedge,Nrows);
else
    xc = linspace(box(1)+boxedge,box(2)-boxedge,Ncols);
    yc = linspace(box(3)+boxedge,box(4)-boxedge,Nrows);
end

% store initial guess locations (centers of ZOI)
[Xc, Yc] = meshgrid(xc,yc);

% create ZOI corner coordinate vectors
zoiX = [Xc(:)-0.5*Lzoi Xc(:)+0.5*Lzoi Xc(:)+0.5*Lzoi Xc(:)-0.5*Lzoi];
zoiY = [Yc(:)-0.5*Lzoi Yc(:)-0.5*Lzoi Yc(:)+0.5*Lzoi Yc(:)+0.5*Lzoi];
searchX = [Xc(:)-0.5*Lsearch Xc(:)+0.5*Lsearch Xc(:)+0.5*Lsearch Xc(:)-0.5*Lsearch];
searchY = [Yc(:)-0.5*Lsearch Yc(:)-0.5*Lsearch Yc(:)+0.5*Lsearch Yc(:)+0.5*Lsearch];

% store for plotting later on
iguess.subset.X = zoiX';
iguess.subset.Y = zoiY';

% for each image
for id = 1:N
    bcwaitbar(H,id/N,sprintf('correlating increment (%d/%d)',id,N));
    
    % for each initial guess location
    for kz = 1:Nzoi
        if id == 1
            % if first image, iguess = 0
            iguess.x(Niguess+kz,id) = Xc(kz);
            iguess.y(Niguess+kz,id) = Yc(kz);
            iguess.ux(Niguess+kz,id) = 0;
            iguess.uy(Niguess+kz,id) = 0;
        else
            % get pixel indices for the search ZOI
            Ima = find(x >= searchX(kz,1) & x <= searchX(kz,2));
            Ina = find(y >= searchY(kz,2) & y <= searchY(kz,3));
            
            % get second image (g)
            xg = x(Ima);
            yg = y(Ina);
            g = Img{id}(Ina,Ima);
            [n, m] = size(g);
            
            % fill f with correct pixels (placed in the center, with
            % boundary of zeros)
            if updated
                % using the previous image as reference
                f = Img{id-1}(y >= zoiY(kz,2) & y <= zoiY(kz,3),x >= zoiX(kz,1) & x <= zoiX(kz,2)); 
            else
                % using the first image as reference
                f = Img{1}(y >= zoiY(kz,2) & y <= zoiY(kz,3),x >= zoiX(kz,1) & x <= zoiX(kz,2)); 
            end
            
            % zero-normalize f
            meanf = mean(f(:));
            stdf = std(f(:));
            ff = (f-meanf)./stdf;

            % zero-normalize g
            meang = mean(g(:));
            stdg = std(g(:));
            g = (g-meang)./stdg;
            
            
            % assign f in the zeropadded large image
            f = zeros(n,m);
            f(yg >= zoiY(kz,2) & yg <= zoiY(kz,3),xg >= zoiX(kz,1) & xg <= zoiX(kz,2)) = ff;
            

            % Cross-Correlation Function (ccf) displacement space, these
            % coordinates are very typical because of the way matlab
            % returns the CCF
            if updated
                ux = 1:0.5:m;
                uy = 1:0.5:n;
            else
                ux = 1:m;
                uy = 1:n;
            end
            ux = ux - ux(floor(length(ux)/2)+1);
            uy = uy - uy(floor(length(uy)/2)+1);
            [Ux, Uy] = meshgrid(ux,uy);
            Ux = ifftshift(Ux);
            Uy = ifftshift(Uy);
            
            % interpolate f and g for the updated case (add one pixel
            % inbetween each pixel). This is to enhance accuracy to
            % mitigate the accumulation of errors in the update method.
            if updated
                f = interp2(f,1,'linear');
                g = interp2(g,1,'linear');
            end
            
            % normalize the image
            scale = sqrt( sum(f(:).^2) * sum(g(:).^2) );

            % Cross-Correlation function
            ccf = ifft2( conj(fft2(f)) .* fft2(g) ) ./ scale;
            
            % find the maximum in the ccf (location of best correlation)
            [cc,I] = max(ccf(:));
            
            % store displacement belonging to the maximum
            if updated
                dux = Ux(I);
                duy = Uy(I);
                
                zoiX(kz,:) = zoiX(kz,:) + dux;
                zoiY(kz,:) = zoiY(kz,:) + duy;
                searchX(kz,:) = searchX(kz,:) + dux;
                searchY(kz,:) = searchY(kz,:) + duy;
                
                iguess.ux(Niguess+kz,id) = iguess.ux(Niguess+kz,id-1) + dux;
                iguess.uy(Niguess+kz,id) = iguess.uy(Niguess+kz,id-1) + duy;
                iguess.cc(Niguess+kz,id) = cc;
            else
                iguess.ux(Niguess+kz,id) = Ux(I);
                iguess.uy(Niguess+kz,id) = Uy(I);
                iguess.cc(Niguess+kz,id) = cc;
            end
        end
    end
    
    % update the application data
    D.iguess = iguess;
    guidata(H,D);

    % update the figure
    set(S.iguessid,'String',num2str(id));
    set(S.iguessslider,'Value',id);
    drawnow
    plotiguess(H)
    
end

% test if new guesses are reasonable, delete an initial guess point if the
% displacement any any frame was bigger than 30% of the search box (30%
% seen from the center means the point moved to within 20% of the edge)
Lreject = D.gui.ini.iguessreject.value*Lsearch;
if updated
    failx = sum(abs(diff(D.iguess.ux(end-Nzoi+1:end,:),[],2)) > Lreject,2);
    faily = sum(abs(diff(D.iguess.uy(end-Nzoi+1:end,:),[],2)) > Lreject,2);
else
    failx = sum(abs(D.iguess.ux(end-Nzoi+1:end,:)) > Lreject,2);
    faily = sum(abs(D.iguess.uy(end-Nzoi+1:end,:)) > Lreject,2);
end
I = find(failx > 0 | faily > 0);
D.iguess.x(Niguess+I,:) = [];
D.iguess.y(Niguess+I,:) = [];
D.iguess.ux(Niguess+I,:) = [];
D.iguess.uy(Niguess+I,:) = [];

% update info
info(1).str = 'Niguess';
info(1).val = length(D.iguess.x);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Auto initial guess routine done');
bcwaitbar(H);

% update the application data
guidata(H,D);

% enable gui controls
set(D.gui.waithandles,'Enable','on');

drawnow
plotiguess(H)


% ==================================================
function [] = iguessaddmanual(varargin)
% Button Callback: Add manual initial guess point
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

if isfield(D.files,'image')
    ImgA = {D.files.image};
else
    return
end
if isfield(D.files,'imageproc')
    ImgB = {D.files.imageproc};
elseif isfield(D.files,'image')
    ImgB = {D.files.image};
else
    return
end

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Manual initial guess routine started');

% get image properties
N = length(ImgA);
[n, m] = size(ImgA{1});
x = 1:m;
y = 1:n;

% get previous guesses
if isfield(D,'iguess')
    iguess = D.iguess;
    Niguess = length(iguess.x);
else
    iguess.x = [];
    iguess.y = [];
    iguess.ux = [];
    iguess.uy = [];
    Niguess = 0;
end

% get zoom options
zoomsize = str2double(get(S.iguesszoomsize,'String'));
zoomselect = get(S.iguesszoomselect,'Value');

% % disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% initialize the reference coordinate
refx = [];
refy = [];

% for each image frame
for k = 1:N
    % update the figure
    set(S.iguessid,'String',num2str(k));
    set(S.iguessslider,'Value',k);
    drawnow
    
    % plot the four images
    for ka = 1:4
        if ka == 1
            A = ImgA{1};
            ha = S.axes41;
            id = 1;
            titlestr = sprintf('id:%d %s',1,'original');
            tag = 'iguess_image1';
        elseif ka == 2
            A = ImgA{k};
            ha = S.axes42;
            id = k;
            titlestr = sprintf('id:%d %s',k,'original');
            tag = 'iguess_image2';
        elseif ka == 3
            A = ImgB{1};
            ha = S.axes43;
            id = 1;
            titlestr = sprintf('id:%d %s',1,'processed');
            tag = 'iguess_image3';
        elseif ka == 4
            A = ImgB{k};
            ha = S.axes44;
            id = k;
            titlestr = sprintf('id:%d %s',k,'processed');
            tag = 'iguess_image4';
        end
        
        % change the image
        set(H,'CurrentAxes',ha)
        hi = findobj(H,'Tag',tag);
        set(hi,'CData',A);
        title(titlestr)
        
        htxt = text(0.01,0.03,sprintf('%d/%d',id,N),'units','normalized');
        set(htxt,'color',D.gui.color.fg)
        set(htxt,'FontSize',D.gui.fontsize(3))
        set(htxt,'FontWeight','bold')
    end
    
    % change the border color, to indicate where to look
    if k == 1
        % first image. look left top
        set(H,'CurrentAxes',S.axes42)
        set(gca,'Xcolor','k')
        set(gca,'Ycolor','k')
        set(H,'CurrentAxes',S.axes41)
        set(gca,'Xcolor',D.gui.color.axeshl)
        set(gca,'Ycolor',D.gui.color.axeshl)
    else
        % consecutive images, look right top
        set(H,'CurrentAxes',S.axes41)
        set(gca,'Xcolor','k')
        set(gca,'Ycolor','k')
        set(H,'CurrentAxes',S.axes42)
        set(gca,'Xcolor',D.gui.color.axeshl)
        set(gca,'Ycolor',D.gui.color.axeshl)
    end
    
    % if using zoom select, zoom the image
    if zoomselect
        if k == 1
            % for the first image, select where to zoom
            ht = title('Select zoom area (i.e. a material point)');
            set(ht,'color',D.gui.color.axeshl);
            
            % using own variant of ginput, to create a more visible cursor
            [xc, yc] = ginputjn(1);
            xlim = xc + [-0.5 0.5]*zoomsize;
            ylim = yc + [-0.5 0.5]*zoomsize;
        else
            % for other images, don't ask, but center the FOV on the
            % previous location
            set(ht,'color',D.gui.color.axeshl);
            xlim = lastx + [-0.5 0.5]*zoomsize;
            ylim = lasty + [-0.5 0.5]*zoomsize;
        end
        
        if k == 1
            % for the first image only zoom the left axes
            set(S.axes41,'xlim',xlim);
            set(S.axes41,'ylim',ylim);
            set(S.axes43,'xlim',xlim);
            set(S.axes43,'ylim',ylim);
        else
            % for the other images zoom the right axes
            set(S.axes42,'xlim',xlim);
            set(S.axes42,'ylim',ylim);
            set(S.axes44,'xlim',xlim);
            set(S.axes44,'ylim',ylim);
            
            % center the zoom around the reference location for the left
            % axes
            xlim = refx + [-0.5 0.5]*zoomsize;
            ylim = refy + [-0.5 0.5]*zoomsize;
            set(S.axes41,'xlim',xlim);
            set(S.axes41,'ylim',ylim);
            set(S.axes43,'xlim',xlim);
            set(S.axes43,'ylim',ylim);
        end
    end

    % prepare axes for material point selection
    ht = title('Select material point');
    set(ht,'color',D.gui.color.axeshl);
    
    % select the material point
    if k == 1
        % if first image, this is the reference
        [refx, refy] = ginputjn(1);
        defx = refx;
        defy = refy;
        defk = k;
        lastx = refx;
        lasty = refy;
        
        % Plot the reference
        href = findobj(H,'Tag','iguess_pointscref');
        if ~isempty(href);
            set(href,'XData',refx,'YData',refy);
        else
            set(H,'CurrentAxes',S.axes41)
            set(gca,'NextPlot','add')
            plot(refx,refy,'*','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscref')
            set(H,'CurrentAxes',S.axes42)
            set(gca,'NextPlot','add')
            plot(refx,refy,'*','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscref')
            set(H,'CurrentAxes',S.axes43)
            set(gca,'NextPlot','add')
            plot(refx,refy,'*','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscref')
            set(H,'CurrentAxes',S.axes44)
            set(gca,'NextPlot','add')
            plot(refx,refy,'*','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscref')
        end
        
        % Plot the deformed
        hdef = findobj(H,'Tag','iguess_pointscdef');
        if ~isempty(hdef);
            set(hdef,'XData',defx,'YData',defy);
        else
            set(H,'CurrentAxes',S.axes42)
            set(gca,'NextPlot','add')
            plot(defx,defy,'.','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscdef')
            set(H,'CurrentAxes',S.axes44)
            set(gca,'NextPlot','add')
            plot(defx,defy,'.','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscdef')
        end
        drawnow
    else
        % if other image, this is a deformed location
        [defx(k), defy(k), but] = ginputjn(1);
        if but == 3
            % if the right button is used, this frame is skipped, store
            % NaNs instead.
            defx(k) = NaN;
            defy(k) = NaN;
            defk(k) = NaN;
        else
            % otherwise store new material point location
            lastx = defx(k);
            lasty = defy(k);
            defk(k) = k;
        end
        
        hdef = findobj(H,'Tag','iguess_pointscdef');
        if ~isempty(hdef);
            set(hdef,'XData',defx,'YData',defy);
        else
            set(H,'CurrentAxes',S.axes42)
            set(gca,'NextPlot','add')
            plot(refx,refy,'.','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscdef')
            set(H,'CurrentAxes',S.axes44)
            set(gca,'NextPlot','add')
            plot(refx,refy,'.','Color',D.gui.color.fg,'Markersize',D.gui.markersize,'Tag','iguess_pointscdef')
            drawnow
        end
        
    end
end

% clear the selection points
delete(findobj(H,'Tag','iguess_pointscref'));
delete(findobj(H,'Tag','iguess_pointscdef'));

% Interpolate skipped frames
% =================

% prepare some data
x = refx;
y = refy;
ux = defx(~isnan(defx)) - x;
uy = defy(~isnan(defy)) - y;
k = defk(~isnan(defk));

% interpolate
if ux < N
    i = 1:N;
    ux = interp1(k,ux,i,'linear','extrap');
    uy = interp1(k,uy,i,'linear','extrap');
end

% reset the axes to the default colors
ha = [S.axes41,S.axes42,S.axes43,S.axes44];
for k = 1:4
    set(get(ha(k),'title'),'Color','k');    
    set(ha(k),'Xcolor','k','Ycolor','k')
    set(ha(k),'Xlim',[1, m],'Ylim',[1, n]);
end

% store the new initial guess points
D.iguess.x(Niguess+1,1) = x;
D.iguess.y(Niguess+1,1) = y;
D.iguess.ux(Niguess+1,:) = ux;
D.iguess.uy(Niguess+1,:) = uy;

% update info
info(1).str = 'Niguess';
info(1).val = length(D.iguess.x);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] Manual initial guess routine done');

% update the application data
guidata(H,D);

plotiguess(H)

% enable gui controls
set(D.gui.waithandles,'Enable','on');

% ==================================================
function [] = iguessframezero(varargin)
% Button Callback: Zero all the points in this frame
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess for all frames,';'create one by pressing <Add auto> or <Add manual>'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end
if isempty(D.iguess.x)
    return
end

% get current frame
id = str2double(get(S.iguessid,'String'));
D.iguess.ux(:,id) = 0;
D.iguess.uy(:,id) = 0;

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[4] Initial guess set to zero for frame %d',id));

% update the application data
guidata(H,D);

plotiguess(H)


% ==================================================
function [] = iguessframeauto(varargin)
% Button Callback: Perform auto iguess for this frame
% this function is very similar to the iguessaddauto, see the comments in
% that function for more info.

H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess for some frames,';'create one by pressing <Add auto> or <Add manual>'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end
if isempty(D.iguess.x)
    return
end
iguess = D.iguess;

% get current frame
id = str2double(get(S.iguessid,'String'));
if id == 1
    return
end

if isfield(D.files,'imageproc')
    Img = {D.files.imageproc};
elseif isfield(D.files,'image')
    Img = {D.files.image};
else
    msgstr = {'This action requires an initial guess for some frames,';'create one by pressing <Add auto> or <Add manual>'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% % disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% coordinate space
[n, m] = size(Img{1});
x = 1:m;
y = 1:n;


% Ask for iguess options
dlg_title = 'Automatic initial guess options';
prompt = {'Subset Size',...
          'Search Size',...
          'Updated'};
def{1} = D.gui.ini.iguesssubsetsize.value;
def{2} = D.gui.ini.iguesssearchsize.value;
def{3} = D.gui.ini.iguessupdated.value;
num_lines = 1;
answer = inputdlgjn(prompt,dlg_title,num_lines,def,dlgposition(H));
if isempty(answer)
    % enable gui controls
    set(D.gui.waithandles,'Enable','on');
    return
end
Lzoi = eval(answer{1});
Lsearch = eval(answer{2});
updated = eval(answer{3});

% store the defaults
D.gui.ini.iguesssubsetsize.value = answer{1};
D.gui.ini.iguesssearchsize.value = answer{2};
D.gui.ini.iguessupdated.value = answer{3};

Xc = D.iguess.x;
Yc = D.iguess.y;
Nzoi = length(Xc);

zoiX = [Xc(:)-0.5*Lzoi Xc(:)+0.5*Lzoi Xc(:)+0.5*Lzoi Xc(:)-0.5*Lzoi];
zoiY = [Yc(:)-0.5*Lzoi Yc(:)-0.5*Lzoi Yc(:)+0.5*Lzoi Yc(:)+0.5*Lzoi];
searchX = [Xc(:)-0.5*Lsearch Xc(:)+0.5*Lsearch Xc(:)+0.5*Lsearch Xc(:)-0.5*Lsearch];
searchY = [Yc(:)-0.5*Lsearch Yc(:)-0.5*Lsearch Yc(:)+0.5*Lsearch Yc(:)+0.5*Lsearch];

for kz = 1:Nzoi
    if mod(kz,5) == 0
        bcwaitbar(H,kz/Nzoi,sprintf('correlating subset (%d/%d)',kz,Nzoi));
    end
    if id == 1
        iguess.x(kz,id) = Xc(kz);
        iguess.y(kz,id) = Yc(kz);
        iguess.ux(kz,id) = 0;
        iguess.uy(kz,id) = 0;
    else
        Ima = find(x >= searchX(kz,1) & x <= searchX(kz,2));
        Ina = find(y >= searchY(kz,2) & y <= searchY(kz,3));
        xg = x(Ima);
        yg = y(Ina);
        g = Img{id}(Ina,Ima);
        [n, m] = size(g);
        f = zeros(n,m);
        
        if updated
            f(yg >= zoiY(kz,2) & yg <= zoiY(kz,3),...
                xg >= zoiX(kz,1) & xg <= zoiX(kz,2))...
                = Img{id-1}(...
                y >= zoiY(kz,2) & y <= zoiY(kz,3),...
                x >= zoiX(kz,1) & x <= zoiX(kz,2));
        else
            f(yg >= zoiY(kz,2) & yg <= zoiY(kz,3),...
                xg >= zoiX(kz,1) & xg <= zoiX(kz,2))...
                = Img{1}(...
                y >= zoiY(kz,2) & y <= zoiY(kz,3),...
                x >= zoiX(kz,1) & x <= zoiX(kz,2));
        end
        
        % ccf displacement space
        if updated
            ux = 1:0.5:m;
            uy = 1:0.5:n;
        else
            ux = 1:m;
            uy = 1:n;
        end
        ux = ux - ux(floor(length(ux)/2)+1);
        uy = uy - uy(floor(length(uy)/2)+1);
        [Ux, Uy] = meshgrid(ux,uy);
        Ux = ifftshift(Ux);
        Uy = ifftshift(Uy);
        
        % interpolate f and g
        if updated
            f = interp2(f,1,'linear');
            g = interp2(g,1,'linear');
        end
        
        mnf = mean(f(:));
        mng = mean(g(:));
        fbar = f - mnf;
        gbar = g - mng;
        scale = sqrt( sum(fbar(:).^2) * sum(gbar(:).^2) );
        
        % Cross-Correlation function
        ccf = ifft2( conj(fft2(fbar)) .* fft2(gbar) ) ./ scale;
        
        % find the maximum in the ccf
        [cc,I] = max(ccf(:));
        
        if updated
            dux = Ux(I);
            duy = Uy(I);
            
            zoiX(kz,:) = zoiX(kz,:) + dux;
            zoiY(kz,:) = zoiY(kz,:) + duy;
            searchX(kz,:) = searchX(kz,:) + dux;
            searchY(kz,:) = searchY(kz,:) + duy;
            
            iguess.ux(kz,id) = iguess.ux(kz,id-1) + dux;
            iguess.uy(kz,id) = iguess.uy(kz,id-1) + duy;
            iguess.cc(kz,id) = cc;
        else
            iguess.ux(kz,id) = Ux(I);
            iguess.uy(kz,id) = Uy(I);
            iguess.cc(kz,id) = cc;
        end
    end
end

% update the application data
D.iguess = iguess;

% test if new guesses are reasonable
Lreject = D.gui.ini.iguessreject.value*Lsearch;
if updated
    failx = sum(abs(diff(D.iguess.ux(end-Nzoi+1:end,:),[],2)) > Lreject,2);
    faily = sum(abs(diff(D.iguess.uy(end-Nzoi+1:end,:),[],2)) > Lreject,2);
else
    failx = sum(abs(D.iguess.ux(end-Nzoi+1:end,:)) > Lreject,2);
    faily = sum(abs(D.iguess.uy(end-Nzoi+1:end,:)) > Lreject,2);
end
I = find(failx > 0 | faily > 0);
D.iguess.x(I,:) = [];
D.iguess.y(I,:) = [];
D.iguess.ux(I,:) = [];
D.iguess.uy(I,:) = [];

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[4] Auto initial guess performed for frame %d',id));
bcwaitbar(H);

% update the application data
guidata(H,D);

plotiguess(H)

% enable gui controls
set(D.gui.waithandles,'Enable','on');


% ==================================================
function [] = iguessframeint(varargin)
% Button Callback: Interpolate the initial guess for this from fom other
% frames
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess for all frames,';'create one by pressing <Add auto> or <Add manual>'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end
if isempty(D.iguess.x)
    return
end
iguess = D.iguess;

% get current frame
id = str2double(get(S.iguessid,'String'));
if id == 1
    return
end

ux = iguess.ux;
uy = iguess.uy;
[N, Nid] = size(ux);

% time vector
ti = 1:Nid;
t  = 1:Nid;

% remove the time values for the current frame
t(id) = [];
% remove the displacements for hte current frame
ux(:,id) = [];
uy(:,id) = [];

% assignin('base','ux',ux)
% assignin('base','uy',ux)
% assignin('base','t',t)
% assignin('base','ti',ti)

% interpolate
for k = 1:N
    if mod(k,5) == 0
        bcwaitbar(H,k/N,sprintf('interpolating subset (%d/%d)',k,N));
    end
    uxi(k,:) = interp1(t,ux(k,:),ti,'linear','extrap');
    uyi(k,:) = interp1(t,uy(k,:),ti,'linear','extrap');
end

D.iguess.ux = uxi;
D.iguess.uy = uyi;

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[4] Initial guess for frame %d interpolated from other frames',id));
bcwaitbar(H);

% update the application data
guidata(H,D);

plotiguess(H)


% ==================================================
function [] = iguessdelone(varargin)
% Button Callback: Delete an initial guess point (nearest one to click)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess,';'create an initial guess in section 4'};
    msgdlgjn(msgstr,dlgposition(H));
    return
end
if isempty(D.iguess.x)
    return
end

iguess = D.iguess;
x = iguess.x;
y = iguess.y;
N = length(x);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% prepare axes
set(H,'CurrentAxes',S.axes41)
ht = title('Select material point');
set(ht,'color',D.gui.color.axeshl);
h = findobj(H,'type','patch');
delete(h);

% select material point
[xc, yc] = ginputjn(1);
set(ht,'color','k');
title('');

% compute the distances from selection to all points
R = sqrt( (x-xc).^2 + (y-yc).^2 );

% find the nearest one
[~, I] = sort(R);
del = I(1);
keep = 1:N;
keep = setdiff(keep,del);

iguess.x  = iguess.x(keep);
iguess.y  = iguess.y(keep);
iguess.ux = iguess.ux(keep,:);
iguess.uy = iguess.uy(keep,:);

D.iguess = iguess;

% update info
info(1).str = 'Niguess';
info(1).val = length(D.iguess.x);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] One initial guess point deleted');

if isempty(iguess.x)
    D = rmfield(D,'iguess');
end

% update the application data
guidata(H,D);

plotiguess(H)

% enable gui controls
set(D.gui.waithandles,'Enable','on');


% ==================================================
function [] = iguessdelbox(varargin)
% Button Callback: Delete initial guess points selected by a box
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess,';'set the initial guess in section 4'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
if isempty(D.iguess.x)
    return
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

iguess = D.iguess;
x = iguess.x;
y = iguess.y;
n = D.files(1).size(1);
m = D.files(1).size(2);
N = length(x);

% Select an area
set(H,'CurrentAxes',S.axes41);
h = title('position the rectangle, confirm with a doubleclick');
set(h,'color',D.gui.color.axeshl);

% position the box
rect = [0.2*m,0.2*n ; 0.4*m,0.4*n];
A = [rect(1,1),rect(1,2);rect(2,1),rect(1,2);rect(2,1),rect(2,2);rect(1,1),rect(2,2);rect(1,1),rect(1,2)];

% load interactive rectangle tool
A = selectarea(A,'Polyline');

% reset the title
set(h,'color','k');
title('');

del = find(inpolygon(x,y,A(:,1),A(:,2)));
keep = 1:N;
keep = setdiff(keep,del);

iguess.x  = iguess.x(keep);
iguess.y  = iguess.y(keep);
iguess.ux = iguess.ux(keep,:);
iguess.uy = iguess.uy(keep,:);

D.iguess = iguess;

% update info
info(1).str = 'Niguess';
info(1).val = length(D.iguess.x);
D.gui.info = appendinfo(D.gui.info,info);

% update status
stat = sprintf('[4] %d initial guess points deleted',length(del));
D.gui.stat = appendstatus(D.gui.stat,stat);

if isempty(iguess.x)
    D = rmfield(D,'iguess');
end

% update the application data
guidata(H,D);

plotiguess(H)

% enable gui controls
set(D.gui.waithandles,'Enable','on');


% ==================================================
function [] = iguessdelall(varargin)
% Button Callback: Delete all initial guess points
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'iguess')
    msgstr = {'This action requires an initial guess,';'set the initial guess in section 4'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% delete the initial guess
D = rmfield(D,'iguess');

% update info
info(1).str = 'Niguess';
info(1).val = 0;
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[4] All initial guess point deleted');

% update the application data
guidata(H,D);
plotiguess(H);


% =========================================================================
% Section 5: Basis
% =========================================================================

% ==================================================
function [] = basissave(varargin)
% Button Callback: save the basis functions
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uiputfile(fullfile(basedir,'basis.mat'), 'save the basis');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;


% save to file
bgdic.basis    = D.basis;
bgdic.version  = D.version;
bgdic.savetype = 'basis';
save('-v7.3',fullfile(pathname,filename),'bgdic')

% update status
D.gui.stat = appendstatus(D.gui.stat,'[5] Basis saved to file');

% update the application data
guidata(H,D);

set(D.gui.waithandles,'Enable','on');drawnow



% ==================================================
function [] = basisload(varargin)
% Button Callback: load the bases
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uigetfile(fullfile(basedir,'basis.mat'), 'load a basis');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

% load from file
load(fullfile(pathname,filename),'bgdic')
if ~exist('bgdic','var') || ~isfield(bgdic,'savetype')
    msgstr = 'incompatible save file';
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end
if ~strcmp(bgdic.savetype,'basis')
    msgstr = sprintf('incorrect save type (%s), try loading it in another section',bgdic.savetype);
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% load the basis functions
D.basis = bgdic.basis;

if isfield(D,'savetype')
    D = rmfield(D,'savetype');
end

% build the plotphi
if isfield(D,'files') && isfield(D,'roi')
    [n, m] = size(D.files(1).image);
        
    % number of pixels (in x-dir) to use to generate plot versions of the basis
    Nplot = 300;
    plotx = linspace(1,m,Nplot);
    ploty = linspace(1,n,round(n*Nplot/m));

    Nb = length(D.basis);
    for k = 1:Nb
        % computing the plot shapes
        D.basis(k).plotx = plotx;
        D.basis(k).ploty = ploty;
        D.basis(k).plotphi = phibuild(plotx,ploty,D.roi,D.basis(k),H,'single');
        D.basis(k).Nphi = size(D.basis(k).plotphi,2);
    end
end
    

% update status
D.gui.stat = appendstatus(D.gui.stat,'[5] Basis loaded from file');

% update the application data
guidata(H,D); 
drawnow;
basisset(H);

set(D.gui.waithandles,'Enable','on');drawnow

plotbasis([],[],H);



% ==================================================
function [] = phiid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'basis')
    return;
end
k = get(S.basislist,'Value');
Nphi = D.basis(k).Nphi;

% get the updated id
id = str2double(get(S.phiid,'String'));
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nphi]);
% update the gui
set(S.phislider,'Value',id);
plotbasis([],[],H);

% ==================================================
function [] = phislider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'basis')
    return;
end
k = get(S.basislist,'Value');
Nphi = D.basis(k).Nphi;

% get the updated id
id = get(S.phislider,'Value');
% rounding and limits
id = round(id);
id = max([id 1]);
id = min([id Nphi]);
% update the gui
set(S.phiid,'String',num2str(id));
plotbasis([],[],H);

% ==================================================
function [] = basislist(varargin)
% Button Callback: Change the currently active basis set
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% get the selected basis
basislst = get(S.basislist,'String');
id       = get(S.basislist,'Value');

if ~isfield(D,'basis') || (length(D.basis) < id) || isempty(D.basis(id).name)
    set(S.basistype,'Value',D.gui.ini.basistype.value);
    set(S.basisname,'String',sprintf('basis %d',id));
    set(S.basisboundary,'String',D.gui.ini.basisboundary.value);
    set(S.basisrows,'String',D.gui.ini.basisrows.value);
    set(S.basiscols,'String',D.gui.ini.basiscols.value);
    set(S.basisorder,'String',D.gui.ini.basisorder.value);
    basistype([],[],H)
    return
else
    set(S.basistype,'Value',D.basis(id).typeval);
    set(S.basisname,'String',D.basis(id).name);
    set(S.basisboundary,'String',D.basis(id).boundratio);
    set(S.basisrows,'String',D.basis(id).Nrows);
    set(S.basiscols,'String',D.basis(id).Ncols);
    set(S.basisorder,'String',D.basis(id).order);
    basistype([],[],H)
end

% update sliders
phiid = str2double(get(S.phiid,'String'));
if ~isfield(D.basis,'plotphi') || isempty(D.basis(id).plotphi)
    Nphi = 1;
else
    Nphi = size(D.basis(id).plotphi,2);
end
phiid = min([phiid Nphi]);

% setting the maximum slider position
slidermax = max([Nphi 2]);
sliderstep = [1/(slidermax-1) 10/(slidermax-1)];

% update sliders (related to basis functions)
set(S.phiid,'String',num2str(phiid))
set(S.phislider,'Max',slidermax)
set(S.phislider,'Sliderstep',sliderstep)
set(S.phislider,'Value',phiid)

plotbasis([],[],H);


% ==================================================
function [] = basisnew(varargin)
% Button Callback: Add a new basis set
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% get the selected basis
basislst = get(S.basislist,'String');
id       = get(S.basislist,'Value');
Nb = length(basislst);

str = sprintf('basis %d',Nb+1);
basislst{Nb+1} = str;
set(S.basislist,'String',basislst);
set(S.basislist,'Value',Nb+1);
D.basis(Nb+1).name = [];

% initiate defaults
set(S.basistype,'Value',D.gui.ini.basistype.value);
set(S.basisname,'String',sprintf('basis %d',Nb+1));
set(S.basisboundary,'String',D.gui.ini.basisboundary.value);
set(S.basisrows,'String',D.gui.ini.basisrows.value);
set(S.basiscols,'String',D.gui.ini.basiscols.value);
set(S.basisorder,'String',D.gui.ini.basisorder.value);

% update the application data
guidata(H,D);
basistype([],[],H)



% ==================================================
function [] = basisdel(varargin)
% Button Callback: Delete a basis set
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% get the selected basis
basislst = get(S.basislist,'String');
id       = get(S.basislist,'Value');
Nb = length(basislst);

I = setdiff(1:Nb,id);
basislst = basislst(I);
D.basis = D.basis(I);

set(S.basislist,'String',basislst);
if id == 1
    set(S.basislist,'Value',id);
else
    set(S.basislist,'Value',id-1);
end


% update the application data
guidata(H,D);

if isempty(D.basis)
    % if all are deleted, create a new one
    Nb = 0;
    str = sprintf('basis %d',Nb+1);
    basislst{Nb+1} = str;
    set(S.basislist,'String',basislst);
    set(S.basislist,'Value',Nb+1);
    D.basis(Nb+1).name = [];

    % initiate defaults
    set(S.basistype,'Value',D.gui.ini.basistype.value);
    set(S.basisname,'String',sprintf('basis %d',id));
    set(S.basisboundary,'String',D.gui.ini.basisboundary.value);
    set(S.basisrows,'String',D.gui.ini.basisrows.value);
    set(S.basiscols,'String',D.gui.ini.basiscols.value);
    set(S.basisorder,'String',D.gui.ini.basisorder.value);
end

basislst = get(S.basislist,'String');
Nb = length(basislst);

% check if the current selection still makes sense
tags{1,1} = 'prepAbasis';
tags{2,1} = 'prepBbasis';
tags{3,1} = 'prepCbasis';
tags{4,1} = 'finalbasis';
tags{5,1} = 'dicrelbasis';
for k = 1:5
    % get the current string
    strs = get(S.(tags{k}),'String');
    val  = get(S.(tags{k}),'Value');
    val = min([val length(strs)]);
    val = max([val 1]);
    str = strs{val};
    % append to the current list
    if k == 5
        basislst = [{'same as U'} ; basislst];
        Nb = Nb + 1;
    end
    
    % find a matching basis
    id = find(strcmp(str,basislst),1);
    if isempty(id)
        val = min([val Nb]);
        val = max([val 1]);
    else
        val = id;
    end
    set(S.(tags{k}),'Value',val)
    set(S.(tags{k}),'String',basislst)
end


% update the application data
guidata(H,D);

basislist([],[],H);


% ==================================================
function [] = basisduplicate(varargin)
% Button Callback: Copy a basis set
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% get the selected basis
basislst = get(S.basislist,'String');
id       = get(S.basislist,'Value');
Nb = length(basislst);

basislst = [basislst ; basislst(id)];
D.basis = [D.basis , D.basis(id)];
set(S.basislist,'String',basislst);
set(S.basislist,'Value',Nb+1);

% update the application data
guidata(H,D);

basislist([],[],H);


% ==================================================
function [] = basistype(varargin)
% Button Callback: Disable/Enable options upon type selection
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% get the basis type
str = get(S.basistype,'String');
val = get(S.basistype,'Value');
type = str{val};

if strcmp(type,'Polynomial')
    set(S.basisboundary,'Enable','Off');
    set(S.basisrows,'Enable','Off');
    set(S.basiscols,'Enable','Off');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'Harmonic')
    set(S.basisboundary,'Enable','Off');
    set(S.basisrows,'Enable','Off');
    set(S.basiscols,'Enable','Off');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'Zernike')
    set(S.basisboundary,'Enable','Off');
    set(S.basisrows,'Enable','Off');
    set(S.basiscols,'Enable','Off');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'B-Spline')
    set(S.basisboundary,'Enable','On');
    set(S.basisrows,'Enable','On');
    set(S.basiscols,'Enable','On');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'FEM Triangle (a)')
    set(S.basisboundary,'Enable','On');
    set(S.basisrows,'Enable','On');
    set(S.basiscols,'Enable','On');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'FEM Triangle (b)')
    set(S.basisboundary,'Enable','On');
    set(S.basisrows,'Enable','On');
    set(S.basiscols,'Enable','On');
    set(S.basisorder,'Enable','On');
elseif strcmp(type,'FEM Quad')
    set(S.basisboundary,'Enable','On');
    set(S.basisrows,'Enable','On');
    set(S.basiscols,'Enable','On');
    set(S.basisorder,'Enable','On');
end

% ==================================================
function [] = basisbuild(varargin)
% Button Callback: Generate the Mesh
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
n = D.files(1).size(1);
m = D.files(1).size(2);
x = 1:m;
y = 1:n;

% Region of interest
if ~isfield(D,'roi')
    msgstr = {'This action requires a defined ROI,';'set the ROI in section 3'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
roi = D.roi;

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the old basis
if isfield(D,'basis')
    basis = D.basis;
else
    basis(1).type = [];
end

% get the selected basis
basislst = get(S.basislist,'String');
id       = get(S.basislist,'Value');

% get the selected type
typelist = get(S.basistype,'String');
typeval  = get(S.basistype,'Value');
type = typelist{typeval};
basis(id).type = type;
basis(id).typeval = typeval;

% get the other options
boundratio = eval(get(S.basisboundary,'String'));
Nrows = round(eval(get(S.basisrows,'String')));
Ncols = round(eval(get(S.basiscols,'String')));
order = round(eval(get(S.basisorder,'String')));

% set minimum values
boundratio = max([boundratio 0.1]);
Nrows = max([Nrows 2]);
Ncols = max([Ncols 2]);
order = max([order 0]);
set(S.basisboundary,'String',num2str(boundratio));
set(S.basisrows,'String',num2str(Nrows));
set(S.basiscols,'String',num2str(Ncols));
set(S.basisorder,'String',num2str(order))

% store to structure
basis(id).boundratio = boundratio;
basis(id).Nrows = Nrows;
basis(id).Ncols = Ncols;
basis(id).order = order;


% Defining the mesh
% ====================
if strcmp(type,'Polynomial')
    basis(id).type = 'polynomial';
    basis(id).order = order;
    
elseif strcmp(type,'Harmonic')
    basis(id).type = 'harmonic';
    basis(id).order = order;

elseif strcmp(type,'Zernike')
    basis(id).type = 'zernike';
    basis(id).order = order;
    
elseif strcmp(type,'B-Spline')
    basis(id).type = 'bspline';
    basis(id).order = order;

    xknot = boundratiovec(roi(1),roi(2),Ncols,boundratio);
    yknot = boundratiovec(roi(3),roi(4),Nrows,boundratio);

    % store for plotting
    basis(id).xknot = xknot;
    basis(id).yknot = yknot;
    basis(id).Nrows = Nrows;
    basis(id).Ncols = Ncols;

elseif strcmp(type,'FEM Triangle (a)')
    basis(id).type = 'FEM-T';

    % limit polynomial order to 2
    order = max([order 1]);
    order = min([order 2]);
    order = round(order);
    set(S.basisorder,'String',num2str(order))
    basis(id).order = order;

    % generate the mesh
    [nodes,conn] = meshgen_T3a(Nrows,Ncols,boundratio,roi);
    if order == 2
        [nodes,conn] = meshgen_T3toT6(nodes,conn,boundratio,roi);
    end
    
    % store the basis
    basis(id).coordinates = nodes;
    basis(id).connectivity = conn;
    basis(id).Nrows = Nrows;
    basis(id).Ncols = Ncols;

elseif strcmp(type,'FEM Triangle (b)')
    basis(id).type = 'FEM-T';

    % limit polynomial order to 2
    order = max([order 1]);
    order = min([order 2]);
    order = round(order);
    set(S.basisorder,'String',num2str(order))
    basis(id).order = order;

    % generate the mesh
    [nodes,conn] = meshgen_T3b(Nrows,Ncols,boundratio,roi);
    if order == 2
        % convert T3 to T6
        [nodes,conn] = meshgen_T3toT6(nodes,conn,boundratio,roi);
    end
    
    % store the basis
    basis(id).coordinates = nodes;
    basis(id).connectivity = conn;    
    basis(id).Nrows = Nrows;
    basis(id).Ncols = Ncols;

elseif strcmp(type,'FEM Quad')
    basis(id).type = 'FEM-Q';

    % limit polynomial order to 2
    order = max([order 1]);
    order = min([order 2]);
    order = round(order);
    set(S.basisorder,'String',num2str(order))
    basis(id).order = order;

    % generate the mesh
    if order == 1
        [nodes,conn] = meshgen_Q4(Nrows,Ncols,boundratio,roi);
    elseif order == 2
        [nodes,conn] = meshgen_Q8(Nrows,Ncols,boundratio,roi);
    end
    
    % store the basis
    basis(id).coordinates = nodes;
    basis(id).connectivity = conn;    
    basis(id).Nrows = Nrows;
    basis(id).Ncols = Ncols;
end

% delete old basis functions
if isfield(basis(id),'plotphi')
    basis(id).plotphi = [];
end

% number of pixels (in x-dir) to use to generate plot versions of the basis
Nplot = D.gui.ini.basisplotcol.value;
basis(id).plotx = linspace(min(x),max(x),Nplot);
basis(id).ploty = linspace(min(y),max(y),round(n*Nplot/m));

% computing the plot shapes
basis(id).plotphi = phibuild(basis(id).plotx,basis(id).ploty,roi,basis(id),H,'single');
basis(id).Nphi = size(basis(id).plotphi,2);

% add the number of dof to the name
basisname = get(S.basisname,'String');
basislst{id} = sprintf('%s (%d)',basisname,basis(id).Nphi);
% update the name to the list
set(S.basislist,'String',basislst);
basis(id).name = basisname;

% also update the basis names of section 6
tags{1,1} = 'prepAbasis';
tags{2,1} = 'prepBbasis';
tags{3,1} = 'prepCbasis';
tags{4,1} = 'finalbasis';
tags{5,1} = 'dicrelbasis';
for k = 1:5
    if k == 5
        basislst = [{'same as U'} ; basislst];
    end
    set(S.(tags{k}),'String',basislst)
end

% update sliders
Nphi = size(basis(id).plotphi,2);
Nphi = max([Nphi 2]);
set(S.phislider,'Value',1)
set(S.phiid,'String','1')
set(S.phislider,'Max',max([Nphi 2]))
set(S.phislider,'Sliderstep',[1/Nphi 10/Nphi])

% make the indicator green
set(S.secind5,'BackgroundColor',D.gui.color.on)

% update info
info(1).str = sprintf('basis %d name',id);
info(1).val = basis(id).name;
info(2).str = sprintf('basis %d type',id);
info(2).val = basis(id).type;
info(3).str = sprintf('basis %d Nphi',id);
info(3).val = basis(id).Nphi;
D.gui.info = appendinfo(D.gui.info,info);

% update status
stat = sprintf('[5] Basis %d generated',id);
D.gui.stat = appendstatus(D.gui.stat,stat);

% update the application data
D.basis = basis;
set(D.gui.waithandles,'Enable','on');drawnow
guidata(H,D);

plotbasis([],[],H);


% =========================================================================
% Section 6: DIC options
% =========================================================================

% ==================================================
function [] = diclagrangeslider(varargin)
% Button Callback: handle Lagrange/Updated Slider changes

H = varargin{3};
S = guihandles(H);

val = get(S.diclagrangeslider,'Value');
val = round(val*100)/100;
str = num2str(val);

set(S.diclagrangeslider,'Value',val);
set(S.diclagrange,'String',str)


% ==================================================
function [] = dicdimensions(varargin)
% Button Callback: If 3D disable relaxation, if 2D enable

H = varargin{3};
S = guihandles(H);

if get(S.dicdimensions,'Value') == 1 %2D
    set(S.dicrelaxation,'Enable','on')
    if get(S.dicrelaxation,'Value') == 1 %none
        set(S.dicrelbasis,'Enable','off')
    else % brightness or brightness+contrast
        set(S.dicrelbasis,'Enable','on')
    end
elseif get(S.dicdimensions,'Value') == 2 %3D
    set(S.dicrelaxation,'Enable','off')
    set(S.dicrelbasis,'Enable','off')
    set(S.resstraindef,'Value',5)
end

% ==================================================
function [] = dicrelaxation(varargin)
% Button Callback: no relaxation disable basis

H = varargin{3};
S = guihandles(H);

if get(S.dicrelaxation,'Value') == 1 %none
    set(S.dicrelbasis,'Enable','off')
else % brightness or brightness+contrast
    set(S.dicrelbasis,'Enable','on')
end


% ==================================================
function [] = cgon(varargin)
% Button Callback: Switch CG step on

H = varargin{3};
D = guidata(H);
S = guihandles(H);

id = varargin{4};

cgstr{1,1} = 'prepA';
cgstr{2,1} = 'prepB';
cgstr{3,1} = 'prepC';
cgstr{4,1} = 'final';

if id <= 3
    set(S.([cgstr{id} 'on']),'BackgroundColor',D.gui.color.hl);
    set(S.([cgstr{id} 'off']),'BackgroundColor',D.gui.color.bg);
    set(S.([cgstr{id} 'on']),'Value',1);
    set(S.([cgstr{id} 'off']),'Value',0);
end

enable = 'on';
set(S.([cgstr{id} 'blur']),'Enable',enable);
set(S.([cgstr{id} 'level']),'Enable',enable);
set(S.([cgstr{id} 'improcess']),'Enable',enable);
set(S.([cgstr{id} 'basis']),'Enable',enable);
set(S.([cgstr{id} 'convcrit']),'Enable',enable);
set(S.([cgstr{id} 'maxit']),'Enable',enable);
set(S.([cgstr{id} 'gradient']),'Enable',enable);
set(S.([cgstr{id} 'tikhpar1']),'Enable',enable);
set(S.([cgstr{id} 'tikhpar2']),'Enable',enable);
set(S.([cgstr{id} 'tikhsteps']),'Enable',enable);

% ==================================================
function [] = cgoff(varargin)
% Button Callback: Switch CG step off

H = varargin{3};
D = guidata(H);
S = guihandles(H);

id = varargin{4};

cgstr{1,1} = 'prepA';
cgstr{2,1} = 'prepB';
cgstr{3,1} = 'prepC';
cgstr{4,1} = 'final';

if id <= 3
    set(S.([cgstr{id} 'off']),'BackgroundColor',D.gui.color.hl);
    set(S.([cgstr{id} 'on']),'BackgroundColor',D.gui.color.bg);
    set(S.([cgstr{id} 'off']),'Value',1);
    set(S.([cgstr{id} 'on']),'Value',0);
end

enable = 'off';
set(S.([cgstr{id} 'blur']),'Enable',enable);
set(S.([cgstr{id} 'level']),'Enable',enable);
set(S.([cgstr{id} 'improcess']),'Enable',enable);
set(S.([cgstr{id} 'basis']),'Enable',enable);
set(S.([cgstr{id} 'convcrit']),'Enable',enable);
set(S.([cgstr{id} 'maxit']),'Enable',enable);
set(S.([cgstr{id} 'gradient']),'Enable',enable);
set(S.([cgstr{id} 'tikhpar1']),'Enable',enable);
set(S.([cgstr{id} 'tikhpar2']),'Enable',enable);
set(S.([cgstr{id} 'tikhsteps']),'Enable',enable);


% ==================================================
function [] = dicsave(varargin)
% Button Callback: save the dic options 
% this function is also called at the start of correlation see
% the corgo function below
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uiputfile(fullfile(basedir,'dicoptions.ini'), 'save the dic options');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;


% gather settings
ini = iniget(H);

% save to file
iniwrite(fullfile(pathname,filename),ini);

set(D.gui.waithandles,'Enable','on');drawnow

% update the status
D.gui.stat = appendstatus(D.gui.stat,'[6] DIC options saved to file');

% update the application data
guidata(H,D);



% ==================================================
function [] = dicload(varargin)
% Button Callback: load dic options
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uigetfile(fullfile(basedir,'dicoptions.ini'), 'load the dic options');
if isequal(filename,0) || isequal(pathname,0)
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

% set section indicator on
set(S.secind6,'BackgroundColor',D.gui.color.on)
 
% update status
D.gui.stat = appendstatus(D.gui.stat,'[6] DIC options saved to file');
 
% update the application data
guidata(H,D);drawnow;
iniset(H);

set(D.gui.waithandles,'Enable','on');drawnow

% =========================================================================
% Section 7: Correlate
% =========================================================================

% ==================================================
function [] = corgo(varargin)
% Button Callback: Start the corrlation

H = varargin{3};

D = guidata(H);
S = guihandles(H);

% remove any old results
if isfield(D,'res')
    D = rmfield(D,'res');
end

set(S.corgo,'BackgroundColor',D.gui.color.hl);
set(S.corstop,'BackgroundColor',D.gui.color.bg);
set(S.corgo,'Value',1);
set(S.corstop,'Value',0);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% enable stop button
set(S.corstop,'Enable','on');drawnow
% enable continue button
set(S.corcontinue,'Enable','on');drawnow
% set indicator color
set(S.secind7,'BackgroundColor',D.gui.color.off)
% set continue button
set(S.corcontinue,'Value',0)

% update the application data
D.gui.ini = iniget(H);
guidata(H,D); drawnow;

corsequence(H); drawnow;

set(D.gui.waithandles,'Enable','on');drawnow
set(S.secind7,'BackgroundColor',D.gui.color.on)

corstop([],[],H)


% ==================================================
function [] = corstop(varargin)
% Button Callback: Stop the correlation (ASAP)

H = varargin{3};
D = guidata(H);
S = guihandles(H);

set(S.corgo,'BackgroundColor',D.gui.color.bg);
set(S.corstop,'BackgroundColor',D.gui.color.hl);
set(S.corgo,'Value',0);
set(S.corstop,'Value',1);


% ==================================================
function [] = corclear(varargin)
% Button Callback: Clear the correlation data for all increments
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if isfield(D,'cor')
    Ninc = length(D.cor);
    for inc = 1:Ninc
        D.cor(inc).done = 0;
        D.cor(inc).U1 = [];
        D.cor(inc).U2 = [];
        D.cor(inc).U3 = [];
        D.cor(inc).U4 = [];
        D.cor(inc).U5 = [];
        D.cor(inc).U6 = [];
    end
end
delete(findobj(H,'Tag','cor_quiver'));

% update status
set(S.corstatus,'String',{''});drawnow

% update status
D.gui.stat = appendstatus(D.gui.stat,'[7] Correlation results cleared for all increments');

% remove any old results
if isfield(D,'res')
    D = rmfield(D,'res');
end

% update the application data
guidata(H,D);

plotcor(H);

% ==================================================
function [] = corclearinc(varargin)
% Button Callback: Clear the correlation data for one increment
% Actually, nothing is deleted, only the done parameter is set to 0
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% get the updated id
inc = str2double(get(S.corid,'String'));

% don't reset the zero increment
if inc == 0
    return
end

if isfield(D,'cor')
    D.cor(inc).done = 0;
    D.cor(inc).U1 = [];
    D.cor(inc).U2 = [];
    D.cor(inc).U3 = [];
    D.cor(inc).U4 = [];
    D.cor(inc).U5 = [];
    D.cor(inc).U6 = [];
end

stat = get(S.corstatus,'String');
stat = [sprintf('increment (%d) cleared',inc) ; stat];
set(S.corstatus,'String',stat);drawnow

% update status
stat = sprintf('[7] Correlation results cleared for increment %d',inc);
D.gui.stat = appendstatus(D.gui.stat,stat);

% remove any old results
if isfield(D,'res')
    D = rmfield(D,'res');
end

% update the application data
guidata(H,D);

plotcor(H);

% ==================================================
function [] = correstart(varargin)
% Button Callback: Clear the correlation data for all increments
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if ~isfield(D,'cor')
    return
end

Ninc = length(D.cor);
for inc = 1:Ninc
    D.cor(inc).done = 0;
end

% update status
set(S.corstatus,'String',{'(Init. option changed to reuse)';'restarting all'});drawnow

% set the Init. options to reuse
set(S.dicinciguess,'Value',5);drawnow

% update status
D.gui.stat = appendstatus(D.gui.stat,'[7] Restart set for all increments');

% remove any old results
if isfield(D,'res')
    D = rmfield(D,'res');
end

% update the application data
guidata(H,D);

plotcor(H);

% ==================================================
function [] = correstartinc(varargin)
% Button Callback: Clear the correlation data for one increment
% Actually, nothing is deleted, only the done parameter is set to 0
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% get the updated id
inc = str2double(get(S.corid,'String'));

% don't reset the zero increment
if inc == 0
    return
end

if isfield(D,'cor')
    D.cor(inc).done = 0;
end

stat = get(S.corstatus,'String');
stat = [sprintf('restarting increment (%d)',inc) ; stat];
stat = ['(Init. option changed to reuse)' ; stat];
set(S.corstatus,'String',stat);drawnow

% update status
stat = sprintf('[7] Restart set for increment %d',inc);
D.gui.stat = appendstatus(D.gui.stat,stat);

% remove any old results
if isfield(D,'res')
    D = rmfield(D,'res');
end

% update the application data
guidata(H,D);

plotcor(H);

% ==================================================
function [] = corid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.corid,'String'));
% rounding and limits
id = round(id);
id = max([id 0]);
id = min([id Nim-1]);
% update the gui
set(S.corslider,'Value',id);
plotcor(H)

% ==================================================
function [] = corslider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.corslider,'Value');
% rounding and limits
id = round(id);
id = max([id 0]);
id = min([id Nim-1]);
% update the gui
set(S.corid,'String',num2str(id));
plotcor(H)


% ==================================================
function [] = cor2iguess(varargin)
% Button Callback: Convert the correlation to initial guess points
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end
if ~isfield(D,'cor')
    msgstr = {'This action requires a succesfull correlation,';'perform a correlation in section 7'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% ask for grid size
prompt = {'number of points in x direction:','number of points in y direction:'};
dlg_title = 'Select initial guess grid size';
num_lines = 1;
def = {D.gui.ini.cor2iguesscols.value,D.gui.ini.cor2iguessrows.value};
answer = inputdlgjn(prompt,dlg_title,num_lines,def,dlgposition(H));
if isempty(answer)
    return
end

% get previous guesses
if isfield(D,'iguess')
    iguess = D.iguess;
    Niguess = length(iguess.x);
else
    iguess.x = [];
    iguess.y = [];
    iguess.ux = [];
    iguess.uy = [];
    Niguess = 0;
end

Nx = round(str2double(answer{1}));
Ny = round(str2double(answer{2}));

Nx = max([1 Nx]);
Ny = max([1 Ny]);
N = Nx * Ny;

% locations of the iguess points
x = linspace(D.cor(1).xroi(3),D.cor(1).xroi(end-2),Nx);
y = linspace(D.cor(1).yroi(3),D.cor(1).yroi(end-2),Ny);
[X, Y] = meshgrid(x,y);

% initiate masked iguess counter
masked = zeros(N,1);

% get maskimage
[n, m] = size(D.cor(1).U1);
maskimg = double(maskimage(D.cor(1).Imask,n,m));

Ncg = 4;
Ninc = length(D.cor);
for inc = 0:Ninc
    inc1 = 1;
    if inc == 0 || isempty(D.cor(inc).done) || (D.cor(inc).done ~= Ncg);
        % zero increment
        iguess.x(Niguess+1:Niguess+N,1) = X(:);
        iguess.y(Niguess+1:Niguess+N,1) = Y(:);
        iguess.ux(Niguess+1:Niguess+N,inc+1) = zeros(N,1);
        iguess.uy(Niguess+1:Niguess+N,inc+1) = zeros(N,1);
    else
        % other increments
        ux = interp2(D.cor(inc1).xroi,D.cor(inc1).yroi',D.cor(inc).U1,X,Y,'linear',0);
        uy = interp2(D.cor(inc1).xroi,D.cor(inc1).yroi',D.cor(inc).U2,X,Y,'linear',0);
        iguess.ux(Niguess+1:Niguess+N,inc+1) = gather(ux(:));
        iguess.uy(Niguess+1:Niguess+N,inc+1) = gather(uy(:));
        
        % interpolate maskimage on iguess locations
        msk = interp2(D.cor(inc1).xroi,D.cor(inc1).yroi',maskimg,X,Y,'linear',1);
        % store mask value for each inc
        masked = masked + msk(:);
    end
    bcwaitbar(H,inc/Ninc,sprintf('initial guess increment (%d/%d)',inc,Ninc));
end

% remove initial guesses which are in the masked area
I = find(masked ~= 0);
iguess.x(Niguess+I,:) = [];
iguess.y(Niguess+I,:) = [];
iguess.ux(Niguess+I,:) = [];
iguess.uy(Niguess+I,:) = [];

% store results
D.iguess = iguess;

% update info
info(1).str = 'Niguess';
info(1).val = length(D.iguess.x);
D.gui.info = appendinfo(D.gui.info,info);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[7] Correlation results copied to initial guess');
bcwaitbar(H);

% update the application data
guidata(H,D);

plotiguess(H);

% =========================================================================
% Section 8: Results
% =========================================================================

% ==================================================
function [] = resid(varargin)
% Callback for the edit box below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = str2double(get(S.resid,'String'));
% rounding and limits
id = round(id);
id = max([id 0]);
id = min([id Nim-1]);
% update the gui
set(S.resslider,'Value',id);
plotres([],[],H)

% ==================================================
function [] = resslider(varargin)
% Callback for the slider below the figures (in figpanel1)
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    return;
end
Nim = length(D.files);

% get the updated id
id = get(S.resslider,'Value');
% rounding and limits
id = round(id);
id = max([id 0]);
id = min([id Nim-1]);
% update the gui
set(S.resid,'String',num2str(id));
plotres([],[],H)


% ==================================================
function resrenderer(varargin)
% Button Callback: switch renderer
H = varargin{3};
S = guihandles(H);
renderer = get(S.resrenderer,'String');
set(H,'Renderer',renderer{get(S.resrenderer,'Value')});

% ==================================================
function [] = ressave(varargin)
% Button Callback: save the result figure
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% user fileselect
[filename, pathname] = uiputfile(fullfile(basedir,'basicgdic_result.png'), 'Save current result figure as');
if isequal(filename,0) || isequal(pathname,0)
    % if cancel
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

% store the basename
D.gui.ini.basedir.value = pathname;

Hpan = S.figpan8;
Hch = findobj(Hpan,'Type','axes');

set(Hpan,'units','pixels');
set(Hpan,'units','normalized');

% Create the figure
% =====================
pos = [50, 50, D.gui.pos.w, D.gui.pos.h];
Hfig = figure('units','pixels',...
    'position',pos,...
    'name','basic gdic figure');

% copy the object to the new figure window
Hch = copyobj(Hch,Hfig);
colormap(D.gui.color.cmapoverlay);

% fix the colorbar width
Ha = findobj(Hch,'Tag','axes8');
Hc = findobj(Hch,'Tag','Colorbar');
set([Ha Hc],'units','pixels');
hapos = get(Ha,'Position');

pos(1) = hapos(1)+hapos(3)+2;
pos(2) = hapos(2);
pos(3) = 20;
pos(4) = hapos(4);
set(Hc,'Position',pos);
set([Ha Hc],'units','normalized');



% set the Render mode
renderer = get(S.resrenderer,'String');
renderer = renderer{get(S.resrenderer,'Value')};
set(Hfig,'Renderer',renderer)

% fix the extention
filename = [regexprep(filename,'.png$','','ignorecase') '.png'];
savepos = get(Hfig,'Position');
% set the paper position to 1 inch per 100 pixels
set(Hfig,'PaperUnits','inches','PaperPosition',savepos.*[0 0 1e-2 1e-2]);
set(Hfig,'PaperSize',savepos(3:4).*[1e-2 1e-2]);

% save png
if strcmpi(renderer,'OpenGL')
    print(Hfig,fullfile(pathname,filename),'-dpng','-r200','-opengl');
else
    print(Hfig,fullfile(pathname,filename),'-dpng','-r200');
end

set(D.gui.waithandles,'Enable','on');drawnow

% update status
stat = sprintf('[8] Results figure saved (%s)',filename);
D.gui.stat = appendstatus(D.gui.stat,stat);

% update the application data
guidata(H,D);



% ==================================================
function [] = rescrosssection(varargin)
% Button Callback: open the cross-section sub-gui
H = varargin{3};
S = guihandles(H);
D = guidata(H);
if ~isfield(D,'files')
    msgstr = {'This action requires loaded images,';'load images in section 1'};
    msgdlgjn(msgstr,dlgposition(H));
    return;
end


% =========================================================================
% Section 9: Info
% =========================================================================


% ==================================================
function [] = inforeset(varargin)
% Button Callback: Reset the entire GUI
H = varargin{3};
S = guihandles(H);
D = guidata(H);

if isfield(D,'files')
    D = rmfield(D,'files');
end
if isfield(D,'roi')
    D = rmfield(D,'roi');
end
if isfield(D,'iguess')
    D = rmfield(D,'iguess');
end
if isfield(D,'basis')
    D = rmfield(D,'basis');
end
if isfield(D,'cor')
    D = rmfield(D,'cor');
end
if isfield(D,'res')
    D = rmfield(D,'res');
end

% initialize status
str = ['Basic gdic tool started: ' datestr(now,'yyyy-mm-dd')];
D.gui.stat = appendstatus({''},str);
D.gui.stat = appendstatus(D.gui.stat,'rule');

% initialize info
D.gui.info = [];
D.gui.info(1).str = 'date';
D.gui.info(1).val = datestr(now,'yyyy-mm-dd');
D.gui.info(2).str = 'time';
D.gui.info(2).val = datestr(now,'HH:MM:SS');

% Load the stored defaults from file
D.gui.ini = iniread;

% update the gui
guidata(H,D);drawnow
iniset(H);
updateinfo(H)

cla(S.axes1)
cla(S.axes2)
cla(S.axes2hist)
cla(S.axes3)
cla(S.axes4)
cla(S.axes41)
cla(S.axes42)
cla(S.axes43)
cla(S.axes44)
cla(S.axes7)
cla(S.axes8)
    
% ==================================================
function [] = infoguidataput(varargin)
% Button Callback: Put the guidata in the workspace (D)
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% update the gui state
D.gui.ini = iniget(H);

% delete the waithandles since they produce large files in matlab 2014b+
if isfield(D.gui,'waithandles')
    D.gui = rmfield(D.gui,'waithandles');
end

assignin('base','D',D);
evalin('base','D');

% update the waithandles
D.gui.waithandles = waithandlesfun(S);

% update status
D.gui.stat = appendstatus(D.gui.stat,'[9] guidata stored to D');
guidata(H,D);drawnow
updateinfo(H)

% ==================================================
function [] = infoguidataget(varargin)
% Button Callback: Get the guidata from the workspace (D)
H = varargin{3};
S = guihandles(H);
if verLessThan('matlab','8.4')
    evalin('base',['guidata(' num2str(H) ',D)']);
else
    evalin('base',['guidata(' num2str(H.Number) ',D)']);
end

% update status
D = guidata(H);

% set the gui elements
if isfield(D,'files')
    set(S.secind1,'BackgroundColor',D.gui.color.on)
    Nfiles = length(D.files);
    filelist = {D.files(:).name}';
    set(S.filelist,'String',filelist);
    set(S.filelist,'Max',Nfiles);
else
    set(S.secind1,'BackgroundColor',D.gui.color.init)
    Nfiles = 2;
end

% update info
info(1).str = 'Nfiles';
info(1).val = Nfiles;
D.gui.info = appendinfo(D.gui.info,info);

% load the gui state
if ~isfield(D.gui,'ini')
    D.gui.ini = iniread;
end
% set all the uicontrols
iniset(H);

% update the waithandles
D.gui.waithandles = waithandlesfun(S);

D.gui.stat = appendstatus(D.gui.stat,'[9] D loaded to guidata');
guidata(H,D);drawnow


updateinfo(H)



% ==================================================
function [] = infosave(varargin)
% Button Callback: save the entire gui (this file can be used as an
% inputfile)
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

if length(varargin) > 3
    filename = varargin{4};
    if strcmp(filename(1),'/') || strcmp(filename(2),':')
        pathname = '';
    else
        pathname = basedir;
    end
else
    [filename, pathname] = uiputfile(fullfile(basedir,'basicgdic.mat'), 'save the GUI state');
    if isequal(filename,0) || isequal(pathname,0)
        set(D.gui.waithandles,'Enable','on');drawnow
        return
    end
end

% store the basename
D.gui.ini.basedir.value = pathname;

% get the gui state
D.gui.ini = iniget(H);

% update status
D.gui.stat = appendstatus(D.gui.stat,sprintf('[9] basicgdic data saved to %s',fullfile(pathname,filename)));

if length(varargin) == 3
    bcwaitbar(H,0.5,'saving basicgdic data');
end

% save to file
bgdic = D;
bgdic.savetype = 'info';

% delete any results (because they are big and can be regenerated from
% cor)
if isfield(bgdic,'res')
    bgdic = rmfield(bgdic,'res');
end

% delete the waithandles since they produce large files in matlab 2014b+
if isfield(bgdic.gui,'waithandles')
    bgdic.gui = rmfield(bgdic.gui,'waithandles');
end
save('-v7.3',fullfile(pathname,filename),'bgdic')

if length(varargin) == 3
    bcwaitbar(H);
end

set(D.gui.waithandles,'Enable','on');drawnow


guidata(H,D);drawnow



% ==================================================
function [] = infoload(varargin)
% Button Callback: Load the an inputfile (initialize the entire GUI)
H = varargin{3};
S = guihandles(H);
D = guidata(H);

% disable gui controls
set(D.gui.waithandles,'Enable','off');drawnow

% get the previous basedir
basedir = D.gui.ini.basedir.value;

% load defaults
if length(varargin) > 3
    % use the argument
    filename = varargin{4};
    if strcmp(filename(1),'/') || strcmp(filename(2),':')
        % absolute path
        pathname = '';
    else
        % relative path
        pathname = basedir;
    end
else
    % ask for a filename
    [filename, pathname] = uigetfile(fullfile(basedir,'basicgdic.mat'), 'load GUI state');
    if isequal(filename,0) || isequal(pathname,0)
        set(D.gui.waithandles,'Enable','on');drawnow
        return
    end
end

% test if the file exists
if ~exist(fullfile(pathname,filename),'file')
    msgstr = 'save file not found';
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end

if length(varargin) == 3
    bcwaitbar(H,0.5,'loading basicgdic data');
end

% load from file
if isempty(whos('-file',fullfile(pathname,filename),'bgdic'))
    % if autosave, the fields are directly in the file
    bgdic = load(fullfile(pathname,filename));
else
    % normal infosave file
    load(fullfile(pathname,filename),'bgdic')
end
if ~isfield(bgdic,'savetype')
    msgstr = 'incompatible save file';
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end
if ~any(strcmp(bgdic.savetype,{'info','autosave'}))
    msgstr = sprintf('incorrect save type (%s), try loading it in another section',bgdic.savetype);
    msgdlgjn(msgstr,dlgposition(H));
    set(D.gui.waithandles,'Enable','on');drawnow
    return
end
D = bgdic;

% store the basename
D.gui.ini.basedir.value = pathname;

if isfield(D,'savetype')
    D = rmfield(D,'savetype');
end

% update the waithandles
D.gui.waithandles = waithandlesfun(S);

% load the gui state
if ~isfield(D.gui,'ini')
    D.gui.ini = iniread;
end

D.gui.stat = appendstatus(D.gui.stat,sprintf('[9] basicgdic data loaded from %s',fullfile(pathname,filename)));
guidata(H,D);
% set all the uicontrols
iniset(H);

if length(varargin) == 3
    bcwaitbar(H);
end


set(D.gui.waithandles,'Enable','on');drawnow

% ==================================================
function [] = updateinfo(H)
% collect information for the infopanel
S = guihandles(H);
D = guidata(H);

% update info
info(1).str = 'date';
info(1).val = datestr(now,'yyyy-mm-dd');
info(2).str = 'time';
info(2).val = datestr(now,'HH:MM:SS');

% Info panel 1: General Info
info = appendinfo(D.gui.info,info);
N = length(info);
for k = 1:N
    rnam{k,1} = info(k).str;
    data{k,1} = info(k).val;
end
cnam = [];
    
set(S.infogeneral,'Data',data,...
    'RowName',rnam,...
    'ColumnName',cnam);

% update the application data
D.gui.info = info;
guidata(H,D);

% Info panel 2: Correlation Info
if isfield(D,'cor')
    % colnames
    cnam{ 1,1} = 'inc';
    cnam{ 2,1} = 'cg';
    cnam{ 3,1} = 'it';
    cnam{ 4,1} = 'Nphi';
    cnam{ 5,1} = 'Npx';
    cnam{ 6,1} = 'r [%]';
    cnam{ 7,1} = 'du [px]';
    cnam{ 8,1} = 'b';
    cnam{ 9,1} = 'dp';
    cnam{10,1} = 'dr [%]';
    cnam{11,1} = 'tikh';
    cnam{12,1} = 'div';
    cnam{13,1} = 'state';
    
    itable = [];
    Ninc = length(D.cor);
    for inc = 1:Ninc;
        if isfield(D.cor(inc),'cg')
            Ncg = length(D.cor(inc).cg);
            for icg = 1:Ncg
                if isfield(D.cor(inc).cg,'itable')
                    itable = [itable ; D.cor(inc).cg(icg).itable ; zeros(1,13)];
                end
            end
        end
    end
    
    rnam = 'numbered';
    cwidth = {40,40,40,40,50,50,60,60,60,60,60,40,40};
    
    cformat{1,1} = 'short';
    cformat{1,2} = 'short';
    cformat{1,3} = 'short';
    cformat{1,4} = 'short';
    cformat{1,5} = 'short';
    cformat{1,6} = 'short';
    cformat{1,7} = 'short';
    cformat{1,8} = 'short';
    cformat{1,9} = 'short';
    cformat{1,10} = 'short';
    cformat{1,11} = 'short';
    cformat{1,12} = 'short';
    cformat{1,13} = 'short';
    
    set(S.infocorrelation,'Data',itable,...
        'RowName',rnam,...
        'ColumnName',cnam,...
        'ColumnWidth',cwidth,...
        'ColumnFormat',cformat);
end

% Info panel 2: Status
set(S.status,'String',D.gui.stat);


% ==============================
function waithandles = waithandlesfun(S)
% A list of handles fo all controls (buttons etc) which will be temporarily
% disabled during computations. This is important since most computations
% depend on the state of the controls, and thus changing them during the
% computations may lead to strange results.


% Comment those lines of the controls you want to remain always active

waithandles = [];
% waithandles = [waithandles, S.secbut9];
waithandles = [waithandles, S.secbut8];
% waithandles = [waithandles, S.secbut7];
% waithandles = [waithandles, S.secbut6];
% waithandles = [waithandles, S.secbut5];
% waithandles = [waithandles, S.secbut4];
% waithandles = [waithandles, S.secbut3];
% waithandles = [waithandles, S.secbut2];
% waithandles = [waithandles, S.secbut1];

waithandles = [waithandles, S.resslider];
waithandles = [waithandles, S.resid];
waithandles = [waithandles, S.resclim2];
waithandles = [waithandles, S.resclim1];
waithandles = [waithandles, S.corslider];
waithandles = [waithandles, S.corid];
waithandles = [waithandles, S.phislider];
waithandles = [waithandles, S.phiid];
waithandles = [waithandles, S.iguessslider];
waithandles = [waithandles, S.iguessid];
waithandles = [waithandles, S.roislider];
waithandles = [waithandles, S.roiid];
waithandles = [waithandles, S.patslider];
waithandles = [waithandles, S.patid];
waithandles = [waithandles, S.imageslider];
waithandles = [waithandles, S.imageid];

waithandles = [waithandles, S.infoload];
waithandles = [waithandles, S.infosave];
waithandles = [waithandles, S.inforeset];
% waithandles = [waithandles, S.infobut3];
% waithandles = [waithandles, S.infobut2];
% waithandles = [waithandles, S.infobut1];

waithandles = [waithandles, S.rescrosssection];
waithandles = [waithandles, S.ressave];
waithandles = [waithandles, S.resstraindef];
waithandles = [waithandles, S.resunit];
waithandles = [waithandles, S.respixelsize];
waithandles = [waithandles, S.resarrowscaleslider];
waithandles = [waithandles, S.resarrowscale];
waithandles = [waithandles, S.resarrows];
waithandles = [waithandles, S.resalphaslider];
waithandles = [waithandles, S.resalpha];
waithandles = [waithandles, S.resoverlay];
waithandles = [waithandles, S.ressofteningslider];
waithandles = [waithandles, S.ressoftening];
waithandles = [waithandles, S.resunderlay];

% waithandles = [waithandles, S.corstatus];
waithandles = [waithandles, S.corliveview];
waithandles = [waithandles, S.cor2iguess];
waithandles = [waithandles, S.corclear];
waithandles = [waithandles, S.corclearinc];
waithandles = [waithandles, S.correstart];
waithandles = [waithandles, S.correstartinc];
% waithandles = [waithandles, S.corstop];
% waithandles = [waithandles, S.corcontinue];
% waithandles = [waithandles, S.corgo];

waithandles = [waithandles, S.dicload];
waithandles = [waithandles, S.dicsave];
waithandles = [waithandles, S.dicusegpu];
waithandles = [waithandles, S.dicprecision];
waithandles = [waithandles, S.dicmemsave];
waithandles = [waithandles, S.dicautosave];
waithandles = [waithandles, S.diclagrangeslider];
waithandles = [waithandles, S.diclagrange];
waithandles = [waithandles, S.dicinciguess];
waithandles = [waithandles, S.dicmaxdiv];
waithandles = [waithandles, S.dicbestit];
waithandles = [waithandles, S.dicconvparam];
waithandles = [waithandles, S.dicrelbasis];
waithandles = [waithandles, S.dicrelaxation];
waithandles = [waithandles, S.dicdimensions];

waithandles = [waithandles, S.basisload];
waithandles = [waithandles, S.basissave];
waithandles = [waithandles, S.basisalphaslider];
waithandles = [waithandles, S.basisalpha];
waithandles = [waithandles, S.basissoftenslider];
waithandles = [waithandles, S.basissoften];
waithandles = [waithandles, S.basisorder];
waithandles = [waithandles, S.basiscols];
waithandles = [waithandles, S.basisrows];
waithandles = [waithandles, S.basisboundary];
waithandles = [waithandles, S.basisname];
waithandles = [waithandles, S.basistype];
waithandles = [waithandles, S.basisduplicate];
waithandles = [waithandles, S.basisdel];
waithandles = [waithandles, S.basisnew];
waithandles = [waithandles, S.basislist];
waithandles = [waithandles, S.basisdefine];

waithandles = [waithandles, S.iguessload];
waithandles = [waithandles, S.iguesssave];
waithandles = [waithandles, S.iguessdelall];
waithandles = [waithandles, S.iguessdelbox];
waithandles = [waithandles, S.iguessdelone];
waithandles = [waithandles, S.iguesszoomsize];
waithandles = [waithandles, S.iguesszoomselect];
waithandles = [waithandles, S.iguessaddmanual];
waithandles = [waithandles, S.iguessaddauto];
waithandles = [waithandles, S.iguessimprocesstype];
waithandles = [waithandles, S.iguessblur];
waithandles = [waithandles, S.iguessimprocess];

waithandles = [waithandles, S.maskload];
waithandles = [waithandles, S.masksave];
waithandles = [waithandles, S.maskclear];
waithandles = [waithandles, S.maskinvert];
waithandles = [waithandles, S.maskintersect];
waithandles = [waithandles, S.maskdel];
waithandles = [waithandles, S.maskadd];
waithandles = [waithandles, S.masktool];
waithandles = [waithandles, S.roiset];

waithandles = [waithandles, S.patshowACF];
waithandles = [waithandles, S.patshowpat];
waithandles = [waithandles, S.patarea];

waithandles = [waithandles, S.filelist];
waithandles = [waithandles, S.filedown];
waithandles = [waithandles, S.fileup];
waithandles = [waithandles, S.filedel];
waithandles = [waithandles, S.fileadd];

waithandles = [waithandles, S.finaltikhsteps];
waithandles = [waithandles, S.finaltikhpar2];
waithandles = [waithandles, S.finaltikhpar1];
waithandles = [waithandles, S.finalgradient];
waithandles = [waithandles, S.finalmaxit];
waithandles = [waithandles, S.finalconvcrit];
waithandles = [waithandles, S.finalbasis];
waithandles = [waithandles, S.finalimprocess];
waithandles = [waithandles, S.finallevel];
waithandles = [waithandles, S.finalblur];

waithandles = [waithandles, S.prepAtikhsteps];
waithandles = [waithandles, S.prepAtikhpar2];
waithandles = [waithandles, S.prepAtikhpar1];
waithandles = [waithandles, S.prepAgradient];
waithandles = [waithandles, S.prepAmaxit];
waithandles = [waithandles, S.prepAconvcrit];
waithandles = [waithandles, S.prepAbasis];
waithandles = [waithandles, S.prepAimprocess];
waithandles = [waithandles, S.prepAlevel];
waithandles = [waithandles, S.prepAblur];

waithandles = [waithandles, S.prepBtikhsteps];
waithandles = [waithandles, S.prepBtikhpar2];
waithandles = [waithandles, S.prepBtikhpar1];
waithandles = [waithandles, S.prepBgradient];
waithandles = [waithandles, S.prepBmaxit];
waithandles = [waithandles, S.prepBconvcrit];
waithandles = [waithandles, S.prepBbasis];
waithandles = [waithandles, S.prepBimprocess];
waithandles = [waithandles, S.prepBlevel];
waithandles = [waithandles, S.prepBblur];

waithandles = [waithandles, S.prepCtikhsteps];
waithandles = [waithandles, S.prepCtikhpar2];
waithandles = [waithandles, S.prepCtikhpar1];
waithandles = [waithandles, S.prepCgradient];
waithandles = [waithandles, S.prepCmaxit];
waithandles = [waithandles, S.prepCconvcrit];
waithandles = [waithandles, S.prepCbasis];
waithandles = [waithandles, S.prepCimprocess];
waithandles = [waithandles, S.prepClevel];
waithandles = [waithandles, S.prepCblur];

function basicgdic_defaults
filename = 'basicgdic_defaults.ini';
if exist(filename,'file')
    return
end

str = {...
'% These are all the GUI settings for the basicgdic tool.'...
''...
'% This file will be automatically generated if it does not exist, i.e. delete it to reset the settings.'...
''...
'% The settings consist of three parts: TAG:FIELD = VALUE, which correspond to the VALUE of the matlab uicontrol FIELD of the object with handle TAG. Lines starting with a "%", or a "[" are ignored and everything beyond a % or [ is ignored as well. The TAG and FIELD parts are case-insensitive, the VALUE part is usually case-sensitive. The order in which settings appear is not important.'...
''...
'[Directories]'...
'% the basedir option allows the specification of a directory where all save and load '...
'% actions originate. This directory can be relative to the current matlab path, or '...
'% absolute. Loading or saving something from a different path will change the basedir '...
'% to the new path. The other paths (inputfile and autosavefile) can be specified '...
'% relatively from the basedir or absolute. Leave these options empty to specify nothing.'...
'basedir:string = '...
'inputfile:string = '...
'autosavefile:string = '...
''...
'[Init. Guess]'...
'iguessblur:string = 3'...
'iguessimprocesstype:value = 5 [1 none, 2 d|xy|, 3 d|x|, 4 d|y|, 5 d(xy), 6 d(x), 7 d(y)]'...
'iguesszoomselect:value = 1 [0 false, 1 true]'...
'iguesszoomsize:string = 150'...
''...
'iguesssubsetsize:string = 80'...
'iguesssearchsize:string = 160'...
'iguesssubsetrows:string = 5'...
'iguesssubsetcols:string = 5'...
'iguessupdated:string = true'...
''...
'[Basis]'...
'% this section only defines the defaults of an undefined basis, use the Load button in this section to load a full set of basis functions'...
'basistype:value = 4 [Polynomial,Harmonic,Zernike,B-Spline,FEM T3a,FEM T3b,FEM Q4]'...
'basisboundary:string = 1.4142'...
'basiscols:string = 1'...
'basisrows:string = 1'...
'basisorder:string = 2'...
''...
'[Basis plot options]'...
'basissoften:string = 1 % between 0 and 1'...
'basisalpha:string = 0.8 % between 0 and 1'...
''...
'[DIC options]'...
'dicdimensions:value = 1 [2D, Quasi 3D]'...
'dicrelaxation:value = 1 [none,constant,linear,quadratic,cubic]'...
'dicrelbasis:value = 1 [same as U, basis 1, etc.] % note that other bases have not been defined yet'...
''...
'dicconvparam:value = 1 [du,b,dp,dr]'...
'dicbestit:value = 1 [r,du,b,dp,dr,last it]'...
'dicmaxdiv:string = 5'...
'dicinciguess:value = 2 [zero,init. guess,prev. inc.,prev. inc., init. guess, reuse]'...
'diclagrange:string = 0 % between 0 and 1'...
''...
'dicautosave:value = 1 [0,1] % i.e. set to 2 to enable autosaving'...
'dicmemsave:value = 2 [0,1,2] % dont get confused, the value indicates the index of the list, a value of 2 means memsave 1'...
'dicprecision:value = 1 [single,double]'...
'dicusegpu:value = 1 [false,true] % remember, specify the index 1->false, 2->true'...
''...
'[Preparation Step A]'...
'prepAon:value = 0 % 0 for false, 1 for true'...
'prepAblur:string = 10'...
'prepAlevel:string = 3'...
'prepAimprocess:value = 1 [none,d|xy|,d|x|,d|y|,d(xy),d(x),d(y)]'...
'prepAbasis:value = 1 % remember, there are no bases defined yet'...
'prepAconvcrit:string = 1e-2'...
'prepAmaxit:string = 10'...
'prepAgradient:value = 1 [auto,gradfg,gradf,gradg]'...
'prepAtikhpar1:string = 0.001'...
'prepAtikhpar2:string = 1e-8'...
'prepAtikhsteps:string = 5'...
''...
'[Preparation Step B]'...
'prepBon:value = 1 % 0 for false, 1 for true'...
'prepBblur:string = 5'...
'prepBlevel:string = 2'...
'prepBimprocess:value = 1 [none,d|xy|,d|x|,d|y|,d(xy),d(x),d(y)]'...
'prepBbasis:value = 2 % remember, there are no bases defined yet'...
'prepBconvcrit:string = 1e-2'...
'prepBmaxit:string = 20'...
'prepBgradient:value = 1 [auto,gradfg,gradf,gradg]'...
'prepBtikhpar1:string = 0.001'...
'prepBtikhpar2:string = 1e-8'...
'prepBtikhsteps:string = 5'...
''...
'[Preparation Step C]'...
'prepCon:value = 1 % 0 for false, 1 for true'...
'prepCblur:string = 1'...
'prepClevel:string = 1'...
'prepCimprocess:value = 1 [none,d|xy|,d|x|,d|y|,d(xy),d(x),d(y)]'...
'prepCbasis:value = 3 % remember, there are no bases defined yet'...
'prepCconvcrit:string = 1e-3'...
'prepCmaxit:string = 50'...
'prepCgradient:value = 1 [auto,gradfg,gradf,gradg]'...
'prepCtikhpar1:string = 0.0001'...
'prepCtikhpar2:string = 1e-9'...
'prepCtikhsteps:string = 10'...
''...
'[Final Step]'...
'finalblur:string = 0'...
'finallevel:string = 0'...
'finalimprocess:value = 1 [none,d|xy|,d|x|,d|y|,d(xy),d(x),d(y)]'...
'finalbasis:value = 4 % remember, there are no bases defined yet'...
'finalconvcrit:string = 1e-4'...
'finalmaxit:string = 50'...
'finalgradient:value = 1 [auto,gradfg,gradf,gradg]'...
'finaltikhpar1:string = 0.0001'...
'finaltikhpar2:string = 1e-9'...
'finaltikhsteps:string = 10'...
''...
'[Correlate]'...
'corliveview:value = 3 [0, 1, 2] % remember the value represents the index'...
'cor2iguessrows:string = 25'...
'cor2iguesscols:string = 25'...
''...
'[Results]'...
'resunderlay:value = 1 [f,g,gtilde,r,none]'...
'ressoftening:string = 1 % between 0 and 1'...
'resoverlay:value = 1 [U1,U2,U3,U4,U5,U6,exx,eyy,exy,eyx,emaj,emin,eeq,r,q,none]'...
'resalpha:string = 0.8 % between 0 and 1'...
'resarrows:value = 1 [U (x,y),strain (xx,yy),strain (min,maj),none]'...
'resarrowscale:string = 1 % between 0 and 10'...
'resrenderer:value = 1 [OpenGL,Zbuffer,Painters]'...
''...
'respixelsize:string = 1'...
'resunit:value = 5 [nm,um,mm,m,km,px,inch,ft]'...
'resstraindef:value = 1 [small,log,Green-Lagrange,Euler-Almansi,membrane,none]'...
''...
''...
'[Advanced]'...
'% these options do not appear in the GUI'...
'iguessboxedge:value = 5 % the distance in px to the box edge, for the auto iguess points'...
'iguessreject:value = 0.3 % iguess rejected if disp. is bigger (rel. to Lsearch)'...
'cordofsupport:value = 0.1 % the minimum relative area of the support of a basis function'...
'corsparseness:value = 0.3 % use sparse matrices if the average relative support is less'...
'corgradcrit:value = 100 % stop updating auto grad when convcrit < corgradcrit*convcrit'...
'cordivtol:value = 1e-7 % consider diverging if: dr > cordivtol (dr is usually < 0)'...
'basisplotcol:value = 300 % number of pixels to use for the shown bases'...
'patnumfacet:value = 7 % number of facets used for evaluating the pattern'...
'cornquiver:value = 15 % number of arrows to show during correlation'...
};

fid = fopen(filename,'w+t');
fprintf(fid,'%s\n',str{:});
fclose(fid);

function str = basicgdic_help

str = {...
'========================================================================'...
' Basic GDIC tool Help'...
'========================================================================'...
''...
'What is the Basic GDIC tool?'...
'--------------------------------------'...
'This program can be used to perform Global Digital Image Correlation (GDIC), with the goal to extract the displacement (field) which occurred during the capturing of the (sequence of) images. The program is named basic because it only touches on the most general applications of Global DIC. The main purpose of this tool is to provide an informative platform, such that, the user can get a direct sense of the impact of the various options on the correlation process. The programming is oriented to maximize freedom and understanding at the cost of memory and computational efficiency. Especially, in the current implementation, the heavy requirements on the memory will require user awareness.'...
''...
'This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.'...
''...
'Global DIC in general'...
'--------------------------------------'...
'Global DIC refers to DIC methods that solve the displacement field of the entire Region Of Interest (ROI), instead of local Zones Of Interest (ZOI). For this reason the displacement field is approximated with a set of basis functions, (e.g. FEM shape functions) with corresponding degrees of freedom (p). The optimal value for the degrees of freedom is then found by iteratively solving the brightness conservation equation using a Newton-Raphson algorithm. When using a mesh in a FEM simulation, more elements (i.e. more Degrees of Freedom DOF) will give a more accurate result (at the cost of computation time). With DIC this is not always true, more DOF will allow the correlation of a more complex displacement field, however, it will also make the procedure more sensitive to noise, and therefore less accurate. The main take-home-message is then, to use enough DOF required to describe the displacement field of your experiment but use as few as possible to limit noise sensitivity (see also [1,2]). Within the framework of Global DIC it is possible to choose any type of basis with which to discretize the displacement field, this GUI provides a few, selecting the optimum basis is challenging since it depends largely on the case at hand. '...
''...
''...
'Quick Help'...
'--------------------------------------'...
'run basicgdic.m or basicgdic.p with matlab, and the graphical user interface (GUI) will appear.'...
''...
'The basic GDIC tool is organized in three panels:'...
'   Panel 1: selection panel (top left)'...
'   Panel 2: control panel (bottom left)'...
'   Panel 3: figure panel (right)'...
' '...
'Panel 1 is constant, and allows the selection of the 9 sections of the tool. Pressing one of the buttons will switch the other two panels to show the controls and figures related to that particular section. The general procedure is then to start in Section 1 (Images) and work through all the sections ending in Section 8 (Results).'...
''...
'- Section 1: Images'...
'   Select image files using <Add files>'...
'   Put them in order using the <Up> and <Down> buttons'...
'   The first image will be the reference (undeformed)'...
'   Use the slider bar (or mouse wheel) to view the images'...
'   '...
'- Section 2: Pattern Eval.'...
'   This section can be used to evaluate the correlation length of the images, however, it is not required to preform correlation, see more details below.'...
'   '...
'- Section 3: ROI and Mask'...
'   Set a Region Of Interest (ROI) by pressing <Set ROI>'...
'   drag/resize the blue box and confirm by double clicking on the blue box'...
'   All material points within the ROI (in image 1) should remain visible in all other images.'...
'   Masking is optional, see details below.'...
'   '...
'- Section 4: Init. Guess'...
'   Press <Image Processing> to create images which are more friendly for the automatic initial guess routine.'...
'   Press <Add auto> to start the automatic initial guess routine'...
'   Drag/Resize the blue box and confirm by double clicking on the blue box'...
'   Wait while the initial guess points are computed'...
'   Use the slider bar (or mouse wheel) to confirm that the results are acceptable'...
'   Use the delete buttons to remove bad initial guess points.'...
'   '...
'- Section 5: Basis'...
'   Select the first basis <Basis 1>'...
'   Select the basis family (Polynomial, T3, etc) and other options'...
'   Select <Define Basis>'...
'   Use the slider bar (or mouse wheel) to evaluate the produced basis functions'...
'   Preferably, create three basis sets increasing in number of degrees of freedom'...
'   '...
'- Section 6: DIC options'...
'   On the left are general DIC options, the defaults are reasonable, '...
'   see below for more details.'...
'   In the right are four panels related to three preparation steps and a fourth and final'...
'   step. The goal of the preparation steps is to find an initial guess with less demanding'...
'   settings (blurred images and few degrees of freedom). Each consecutive step starts with'...
'   the result of the previous step as initial guess.'...
'   For each step:'...
'      Select the Coarsegrain level N, which will superpixel the image (average N^2 pixels)'...
'      Select the basis used in this CG step in order of complexity,'...
'   '...
'- Section 7: Correlate'...
'   Press <Go> to start the correlation process'...
'   On the right the residual field is shown, this should decrease to zero if no'...
'   acquisition noise is present'...
'   The shown arrows display the displacement field'...
'   Press <Stop> to abort the correlation process'...
'   Pressing <Go> again will continue the previous correlation and skip any increments'...
'   which are already done'...
'   To start fresh, select <Clear All>'...
'   To start from the current results, select <Restart All>'...
'   '...
'- Section 8: Results'...
'   Select a strain definition and wait until the strain fields are computed'...
'   Select a pattern to show (lowest plot level)'...
'   Select an overlay, e.g. x-displacement or strain-xx (intermediate plot level)'...
'   Select a vector field (highest plot level)'...
'   Save to .png using <Save>'...
''...
'- Section 9: Info'...
'   Go back to this help, or show general info or status info produced during the various steps.'...
''...
''...
''...
'========================================================================'...
' Detailed Explanation per Section'...
'========================================================================'...
''...
'Starting the basicgdic application'...
'--------------------------------------'...
'Depending on whether you have the open source version or protected version you will have a basicgdic.m (source) or basicgdic.p (protected) file. If you have the protected version and would like the open version send me an email (see contact information below). The Graphical User Interface will start after running either the basicgdic.m file or the basicgdic.p file with matlab. For ease of use there is a runme.m file delivered with the protected version such that you can double-click it in your file-browser to open it in matlab and then press F5 to run it.'...
''...
'It is also possible to start the application with an inputfile by running basicgdic(''myinputfile.mat''), where myinputfile.mat should be a file such as saved with the <Save> button in the info section (section 9). This will start the GUI and load the inputfile in the same way as it would be loaded by the <Load> button in the same section. Instead of an inputfile the tool also accepts the D structure as input e.g. basisgdic(D).'...
''...
'Optionally, it is possible to run the basicgdic application in headless mode by also specifying an outputfile e.g. basicgdic(''myinputfile.mat'',''myoutputfile.mat''). This will start the GUI in invisible mode, load the inputfile like before, start the correlation process, save everything to myinputfile.mat, and close the application. The file created can be loaded in the GUI by using the <Load> button of section 9 to evaluate the results. The headless mode is intended such, that you can prepare a correlation using the GUI, providing all the necessities such as initial guess and DIC options on a lightweight machine, say your laptop. Then save everything and start the correlation on a more powerful machine, say a cluster, which may not have a display. Currently, the headless mode is not well tested and care should be taken when running on clusters. For instance, matlab will use more CPU cores then strictly claimed, ask for help if you don''t know what this means. Additionally, the basicgdic application will still generate pop-up error messages when encountered, which may cause problems on the cluster, please test on your own machine first.'...
''...
'The basisgdic tool will also return the D structure. However, in normal GUI mode, the output is returned early in the process, i.e. directly after the tool is loaded. Consequently, the contents of D are not very interesting. However, when in headless mode a full D structure is returned, as if the <Guidata to D> button is pressed just before the tool closes.'...
''...
''...
'Changing the basicgdic defaults'...
'--------------------------------------'...
'The basisgdic tool will look for a file named "defaults.ini", stored in the same folder as basicgdic.m. This ini file allows the specification of default settings that will precede the build-in defaults. Be careful, the listed settings must exist in the GUI, i.e. typos will generate an error. Have a look in the lib_bgdic/basicgdic_defaults.ini file for examples, but do not edit that file directly. '...
''...
'The settings consist of three parts: TAG:FIELD = VALUE, which correspond to the VALUE of the matlab uicontrol FIELD of the object with handle TAG. Lines starting with a "%", or a "[" are ignored and everything beyond a % or [ is ignored as well. The TAG and FIELD parts are case-insensitive, the VALUE part is usually case-sensitive. The order in which settings appear is not important.'...
''...
''...
'Section 1. Images'...
'--------------------------------------'...
'The first (or top) image will always be considered as the reference configuration. And the order of appearance will also be considered chronological. Therefore, it is important to sort the files correctly. After loading, the images are shown in the right panel, with the image number shown in the lower left corner. It is possible to view all images after loading by dragging the slider bar (or using the mouse wheel) below the figure. The x and y axes will be in pixels throughout the tool, however in the results section they can be converted to physical dimensions.'...
''...
'It is possible to re-order the images at any stage by revisiting this section. However, if the reference image is changed, by either removing it or moving it to another frame, then all computed results are reset. If new images are added, they are added to the bottom of the list, after which they can be moved to the right location in the list. It should be noted that, new images start with an initial guess of zero. In the initial guess section (section 4) this can be fixed using the <Single frame operations>, see below.'...
''...
'Currently a limited number of file-types are supported (basically those supported by the imread command). After loading, all images are stored as double precision floating point matrices. Extending the tool to allow additional file formats should be relatively easy, have a look at the fileadd function in the basicgdic.m file.'...
''...
''...
'Section 2. Pattern Eval.'...
'--------------------------------------'...
'Besides the image histogram, shown in bottom of the control panel, this section allows the computation of the Auto Correlation Function (ACF) and the correlation length, or correlation radius (zeta). Press <Evaluate> to start the computation, after which a 7x7 grid of circles is shown in the <Show Pattern> window, and the ACF is shown in the <Show ACF> window. The correlation radius is defined as the radial distance at which the ACF intersects 0.5. The correlation radius can be seen as the allowable distance of your initial guess from the real solution in order to have proper convergence. Since the correlation radius depends on the evaluated area, it is different for each degree of freedom, since they all have a different support. A small correlation radius is detrimental for the robustness. However, a small correlation radius also means that the image gradients are large, which results in improved accuracy. This is the reason why DIC methods (and this tool) usually apply a step-wise procedure where first a correlation is performed on a blurred image (or coarse grained image), the result of which is used as the initial guess for a finer image etc. until finally the image is evaluated in full detail.'...
''...
''...
'Section 3. ROI and Mask'...
'--------------------------------------'...
'The Region Of Interest (ROI) is required and defines the area (rectangle) over which the displacement field is computed. Any material point within the ROI in the first image should remain visible in all other images. Otherwise it is not possible to compute a residual (i.e. brightness conservation).'...
''...
'Areas within the ROI which negatively influence the DIC procedure can be excluded from the minimization algorithm by masking them. These are usually areas where there is no sample, no pattern, or a loss of brightness conservation (shadows, reflections, paint loss, etc.).'...
''...
'The displacement field is computed for masked pixels, they just do not have any influence on the solution. If certain degrees of freedom have become ill-supported because most of their supported pixels are masked (more then 80%), then these DOF are also removed from the minimization and the update in those DOF (dp) is set to zero.'...
''...
''...
''...
'Section 4. Init. Guess'...
'--------------------------------------'...
'Most DIC procedures are challenged when large displacements are present. In the internal algorithms the pattern is linearized, and due to the high non-linearity of the pattern, this linearization only has a limited range. Therefore it is important to initialize the DIC routine with a good initial guess. This section allows for generating automatic initial guesses using a local DIC (cross-correlation) algorithm. And if this does not work, a manual selection of material points is also possible.'...
''...
'The first options are related to <Image Processing>, where the images can be blurred using a Gaussian kernel with a chosen standard deviation (sigma), remember that the radius of the kernel is approximately 3*sigma. If the blur value is smaller than 0.3, the image is not blurred. The second option is "coarsegrain", which also blurs the image but now by averaging groups of pixels (N^2 x N^2) into superpixels. This option also reduces the image matrix size, thereby reducing the computational cost. The third option <Im. process> creates a new image by computing the image gradients and combining them.'...
'if f is the image:'...
'none => results in the normal image, i.e. no processing'...
'grad|xy| => abs(df/dx) + abs(df/dy)'...
'grad|x| => abs(df/dx)'...
'grad|y| => abs(df/dy)'...
'grad(xy) => (df/dx) + (df/dy)'...
'grad(x) => (df/dx)'...
'grad(y) => (df/dy)'...
''...
'Pressing <Add auto> starts the automated initial guess routine. As said above a local DIC algorithm is used. After pressing the button a popup dialog appears where the parameters of the method can be set (detailed below). Thereafter a blue rectangle is presented in the top left figure which can be resized and repositioned. The location is confirmed by double clicking on the rectangle. After confirmation, the initial guess points are set and the displacement for each subset is computed for all images.'...
''...
'Add auto parameters:'...
'<Subset size> => the width (and height) of the square subset in the reference image which is searched for in the deformed images.'...
'<Search size> => the width (and height) of the square search area in the deformed images where the reference subset is compared to. The location of the search area is centered around the reference subset location. The search area should be larger than the subset plus the expected displacement. Computed displacements which are at the edge of the search window are automatically discarded.'...
'<Subset rows> => Number of subsets rows'...
'<Subset cols> => Number of subset columns'...
'<Updated> => type "true" to switch the local DIC algorithm between correlating between consecutive images (i.e. updated Lagrange) or "false" to correlate each image with the first image (i.e. total Lagrange).'...
''...
'Afterwards the slider bar (or mouse wheel) can be used to evaluate the generated initial guess points. Any bad initial guess points should be removed before continuing, use the Delete <One>. <Area> or <All> buttons in the bottom of the control panel to delete initial guess points.'...
''...
'If the auto initial guess procedure does not work, a manual procedure is provided by the <Add manual> button. After pressing this button a material point should be selected with the mouse in the top left figure. After which the same material point should be selected for each consecutive image in the top right figure. However some images can be skipped by using the right mouse button. The initial guess will be linearly interpolated/extrapolated for the skipped images.'...
''...
'The box below the <Add manual> button specifies the zoomed area in pixels. If zoomselect is checked (advised) then the first mouse click (in the top left figure) zooms the figures, the second click (in the top left figure) then defines the reference material point. Subsequently, the same material point should be selected in the top right window for each consecutive image as described in the previous paragraph.'...
''...
'The <Single frame operations> allow, as the name suggests, operations in a single frame where:'...
'<Zero> sets the (displacements) initial guess points to zero for the current frame.'...
'<Auto> performs the same automated initial guess routine as above, except only on the current frame.'...
'<Interpolate> interpolates the initial guesses for each point for the current frame based on the initial guesses of all other frames. A linear interpolation scheme in time is used per point, which can also extrapolate. This method assumes a constant time spacing of the image frames, and actually does not use any information from the image itself.'...
''...
''...
'Section 5. Basis'...
'--------------------------------------'...
'The power of Global DIC is in the wide range of choices for basis functions. As discussed before, too few DOF will restrict the kinematics too much resulting in reduced accuracy, and too many DOF will make the method noise sensitive, resulting in reduced accuracy. Consequently there is an optimum, which is the minimum number of DOF which include the (yet unknown) displacement field. Obviously, this minimum depends in the experiment at hand. Less obviously, this minimum also depends on the type of basis, i.e. the family of shape functions.'...
''...
'Usually, a new basis set is defined per coarse grain step, however it is possible to define as many as you like in this section. In most cases 3 coarse grain steps, and thus 3 basis sets, will suffice for the correlation process. All coarse grain steps except the last are only used to as an initial guess for the next coarse grain step. Therefore, their displacement fields (and thus their basis) only has to be good enough to give a adequate initial guess. Consequently, using Polynomial or Harmonic basis sets for the these CG steps tends to work better since these functions are more robust (due to their wide support).'...
''...
'To define a basis set use the <New> button to start a new set, or the <Del> button to remove a set, or make a copy of a set using the <Dupl.> button. Depending on the selected <Type> (more details below) the <Boundary ratio>, <Rows>, <Cols> and <Order> options will be available. The <Basis name> option can be used to give the basis set a custom name. After setting all of these options, the <Define Basis> button can be used to generate the basis functions, which will appear on the right for evaluation. Obviously, it is possible to redefine any of the defined basis sets by just going back to the basis set, changing some options and pressing <Define Basis> again.'...
''...
'Polynomial:'...
'These functions have support over the entire ROI making them highly robust, but not well suited for highly localized deformations. The zero order functions allow constant displacements (i.e. rigid body translation), first order functions allow constant strain, etc. The 2D shapes are computed from the dyadic product of two 1D Legendre polynomials. The polynomials are computed on normalized coordinates which span from -1 to 1 over the ROI in both directions. Consequently, these basis functions are always 1, or -1 in the corners of the ROI.'...
''...
'Harmonic:'...
'Similar to the polynomials, these also have wide support over the entire ROI. They are created by the dyadic product of 1D sine and cosine functions. The zero order is a (non-harmonic) constant function, i.e. a zero order polynomial. Each consecutive order introduces a pair of harmonic functions with one wavelength (a sine and a cosine), with decreasing wavelength for increasing order. The wave lengths are chosen such that they fit "k" times in 1.5 times the width and height of the ROI, centered in the ROI.'...
''...
'Zernike:'...
'These 2D orthogonal pseudo Zernike polynomials are defined on circular coordinates. They are particularly useful for describing lens aberrations.'...
''...
'B-Spline:'...
'With the <Rows> and <Cols> a rectangular grid (control net) is defined in which regular 1D B-Spline functions are used to compute 2D shapes. Zero order B-Splines result in local DIC like behavior, i.e. subsets. First order B-Splines are equal to bi-linear (Q4) FEM shape-functions. And higher order B-Splines are also possible. Strain is the first derivative of the displacement, therefore, to have continues smooth strain fields, second order B-Splines are recommended. The <boundary ratio> option changes the spacing of the boundary knots relative to the spacing of the internal knots. If unsure about which basis fits the best for your case, try the second order B-splines, and play with the number of rows and cols.'...
''...
'FEM Triangle:'...
'This generates a FE mesh of triangular elements. The mesh is generated by initiating a set of nodes in a grid which are connected by Delauney triangulation. The difference between "Triangle (a)" and "Triangle (b)" is the way the nodes are positioned before triangulation. Note that even if the nodes are in a constant grid, that there is some variation in the size of the support for each DOF, depending on the connectivity. The <boundary ratio> option changes the spacing of the boundary nodes relative to the spacing of the internal nodes. This mesh type accepts <order = 1> for linear 3-noded (T3) elements, or <order = 2> for quadratic 6-noded (T6) elements.'...
''...
'FEM Quad:'...
'This generates a FE mesh of quadrilateral elements in a regular grid. The <boundary ratio> option changes the spacing of the boundary nodes relative to the spacing of the internal nodes. This mesh type accepts <order = 1> for bi-linear 4-noded (Q4) elements, or <order = 2> for serendipity 8-noded (Q8) elements. The regular mesh generation is only a limitation of the mesh generator, the GUI can deal with irregular meshes. Although, obtaining the FEM shape-functions requires a reverse mapping operation which can be unstable at the element corners if the elements are too distorted (e.g. more than 15% distortion).'...
''...
'Custom meshes:'...
'The basis functions can be saved to a .mat file, and loaded back into to tool. This also provides a means to load custom meshes. Have a look at the way the .mat file is structured by saving the mesh and loading it in matlab using the "load" command. '...
'For example (assuming a FEM Q4 mesh loaded in basis 1):'...
'   use the <save> button (in section 5) to save the basis.mat file then run:'...
'   >> load basis.mat'...
'   >> bgdic.basis(1).coordinates = bgdic.basis(1).coordinates + 3*randn(size(bgdic.basis(1).coordinates,1),2)'...
'   >> save(''basis.mat'',''bgdic'')'...
'   then load the .mat file using the <load> button (in section 5).'...
''...
'Section 6. DIC Options (General Options)'...
'--------------------------------------'...
'On the left are the general DIC options:'...
''...
'<Dimen.> "2D" or "Quasi 3D", the latter switches to Quasi 3D, which means that the image gray values are interpreted as height values. Use this method if the images are surface topographies.'...
''...
'<Relax.> This option allows the relaxation of the brightness equation, select "none" to disable relaxation. Otherwise, set the desired order of relaxation. The brightness conservation equation will be enriched with a number of terms (equal to the chosen order) i.e. r = f - (g + U3*g^0 + U4*g^1 + U5*g^2 + U6*g^3). It is important to realize that U3,U4,U5 and U6 are fields which require basis functions and degrees of freedom, similarly as the in-plane displacement fields (U1,U2). The basis for the relaxation fields can be set in the next option <Rel. Basis>. Brightness relaxation is essential when there have been irregularities in the image intensity during the experiment. Examples are, changes in lighting conditions or reflections on the sample surface. Typically, "Constant" or "Linear" brightness relaxation is sufficient. '...
''...
'<Rel. Basis> Select the basis set, used to describe the additional brightness (and/or contrast) fields. The "same as U" option uses the same basis set as used for the displacement field, which may vary for each coarse grain step. The other items in the list are basis sets defined in section 5, selecting one of those will cause the same basis set to be used for all coarse grain steps.'...
''...
'<Conv. Cr.> Which parameter to use to test convergence on. "Change in disp." is the maximum of the iterative change in the displacement fields, in pixels. "The right hand member" is the norm of b, which is the object the DIC algorithm is actually reducing to zero. "update in dof" is the norm of dp, which is the iterative change in the degrees of freedom. "change in residual" is the change in the norm of r, in GV. Due to acquisition noise and interpolation errors, the residual usually does not really go to zero, therefore the change in the residual is offered as a convergence parameter.'...
''...
'<Best iter.> After convergence is reached, instead of using the last iteration, it is also possible to use a previous iteration as final result. Select which <Best iter.> parameter to use to identify the best iteration. "residual" uses the iteration with the lowest (mean, absolute) residual. The following four are the same ones as discussed for the convergence criterion. The last option is "last it." simply uses the last computed iteration.'...
''...
'<Max div.> If the correlation procedure causes the residual to increase N consecutive times the coarse grain step is stopped. For each iteration which decreases the residual the counter is reduced by 0.5.'...
''...
'<Next inc.> This option selects how the initial guess is computed for the next increment. "Zero" causes all DOF to be initiated at zero for each increment (only for the first coarse grain step). "init. guess" uses the initial guess points as defined in section 4 for each increment. "prev. inc." uses the displacement field of the previous increment as initial guess for this increment (zero for the first increment). "prev. inc. and init. guess" averages between the second and third option. The "reuse" option re-uses the data currently present in the GUI, i.e. the data from a previous correlation result. This option reverts to "prev. inc." if no previous data for an increment is found.'...
''...
'<Total/Updated Lagrange> A Total DIC scheme always correlates between the first image and the current (i) image. An Updated scheme always correlates between two consecutive images (i-1 and i). The advantage of a total scheme is that it is more accurate since each increment is minimized individually. However, if the pattern changes during the experiment due to shadows, reflection, paint loss, etc. then the current image may be too different from the first image for successful correlation. Consequently an updated scheme is more robust, but errors in the displacement scheme are accumulated over the increments. Beware, using a non-zero value here will make the correlation of the current increment depend on the quality of the correlation of the previous increment. If the previous increment is not converged, the algorithm will skip the current increment, and probably all following ones.'...
''...
'<Auto save> This option causes the correlation results to be saved at the end of each increment. Upon pressing <Go> in section 7, a file selection popup appears to select the location for the autosave file. The file is overwritten after each increment.'...
''...
'<Mem. save> This particular DIC tool is very memory hungry. This is because internally, a matrix L is required to compute M (the correlation matrix) and b (the right hand member). The number of rows in this L matrix is equal to the number of unmasked pixels within the ROI, the number of columns in this L matrix is equal to the number of DOF. For example a 1 Megapixel image with a T3 mesh with 100 nodes will require an L matrix of 1e6*100*2*8 = 1.6 Gb for 2D and even more for Quasi 3D or brightness and contrast relaxation. Choosing <Mem. save = 0> causes exactly this behavior, which is the fastest option if your computer has enough memory. Choosing <Mem. save = 1> partitions the L (and M) matrix into parts separated per dimension (2 for 2D, 3 for 3D or 2D+B, 4 for 2D+B+C). This is slower but halves (or more) the required memory space. Choosing <Mem. save = 2> partitions the L (and M) matrix per DOF, this greatly reduces the memory requirements but these L partitions need to computed over and over since they cannot be stored in memory, this last option is very slow. Instead of using the last option, in most cases using a coarse-grained final image can give faster but slightly les accurate results. Whenever possible (e.g. FEM or B-Spline basis) L is stored sparsely, which will save a lot of memory, however this complicates the prediction of the required memory.'...
''...
'<Precision> This option sets the internal precision to "single" (32-bit) or "double" (64-bit). Setting the precision to single causes the memory requirement to be halved, and greatly improves the computation speed on graphics cards. However, the floating point precision will be reduced from approximately 16 digits to 6 digits behind the decimal. Moreover, Matlab cannot handle single precision sparse matrices, therefore, sparse matrices are disabled when the precision is set to single.'...
''...
'<use GPU> Enabling this transfers some computations to the GPU (Graphics Card) instead of the CPU. The code is not optimized for GPU computing, it merely utilizes the Matlab built-in CUDA capabilities, and therefore requires a CUDA (NVidia) video card. Typically, the video card has less memory than the CPU, which will be limiting when correlating larger images (e.g. 1 Mpx) with many DOF. Additionally, current Matlab versions (2013b) cannot deal with sparse matrices on the GPU, making the problem worse.'...
''...
''...
'Section 6. DIC Options (Preparation Steps and Final Steps)'...
'--------------------------------------'...
'On the right are the options for each step. Three optional preparation steps and one mandatory final step can be configured. The result of each step will serve as the initial guess to the next step. In most DIC cases 2 coarse preparation steps should be sufficient, therefore, the fest step is disabled by default. '...
''...
'The goal of these preparation steps is to perform DIC using less accurate but more forgiving settings. This way the algorithm can converge more easily. The result only needs to be accurate enough to serve as a initial guess for the next step which will have more demanding settings. This mostly boils down to two things: Smoothing the image (blur/coarsegrain) and reducing the number of DOF of the basis (selecting an easier basis e.g. bigger elements).'...
''...
'Blurring an image results in a larger correlation radius, i.e. allows the algorithm to converge without getting stuck in a local minimum. However, blurring also removes information from the image, which reduces the accuracy, but also limits the number of DOF that can be found. The better the provided initial guess, the less coarse graining is required. '...
''...
'Each step has the same options:'...
''...
'<On/Off> Enable, disable this step. '...
''...
'<Blur> Blur the image with a Gaussian kernel with standard radius (sigma) N, remember that the support of the kernel is approximately 3*N. Blurring removes fine detail to allow for longer distance correlation but does not reduce the number of pixels. Consequently, this method has no speed gain, but is less likely to produce singularity problems.'...
''...
'<Coarsegrain> This option creates a new image with superpixels. The superpixels are computed by averaging groups of (N^2 x N^2) pixels. This has the combined effect of blurring and reducing the matrix size. This is the preferred method of coarse graining. Remember that a CG level of 4 will create superpixels of size 16x16 which reduces a 1 Megapixel image to roughly 4000 pixels (62x62), setting the CG level higher is usually not useful.'...
''...
'<Image Process> This option allows the computation of a gradient image, which is then used instead of the real image to correlate with. Gradient images highlight edges, which may be useful for cases with large variations in grayvalues from image to image, e.g. shadows or reflections. Such detrimental features usually have less influence on the gradient images. For most cases this is not recommended.'...
''...
'<Basis> Select one of the bases defined in section 5. Make sure to set the basis with the least DOF for the first step and the basis with the most DOF for the final step.'...
''...
'<Convcrit.> The convergence criteria, i.e. when to consider the step converged. For the first two steps a relatively coarse criteria can be used (e.g. 1e-2) since the obtained results are only initial guesses for the final step.'...
''...
'<Maxit.> Maximum number of iterations in this step. Similar to the <Convcrit.> the first steps do not require perfect convergence, the result only needs to be "good enough" as an initial guess for the next step.'...
''...
'<Image Gradient> This is the image gradient used in the iterative (Newton-Raphson) algorithm, see [3] for more details.'...
'gradf: the classical DIC gradient, the cheapest to compute, but limited to small displacements.'...
'gradg: the consistent DIC gradient, expensive to compute, compatible with large displacements.'...
'gradfg: the average of the above two, expensive to compute, has an advantage when far from the solution, i.e. for bad initial guesses, suited for intermediate displacements (e.g. rotations up to 45°).'...
'auto: starts the routine with gradfg, but stops updating the gradient when close to convergence to reduce the computational cost. The transition occurs when the chosen convergence parameter is smaller than 100 times the set convergence criteria.'...
''...
'<Tikhonov Regularization> This method adds a second potential to the minimization potential which increases when the DOFs move away from their initial guess. The strength of this potential is controlled with the Tikhonov parameter (alpha). It is implemented such that alpha can change in each iteration moving from <alpha(1)> to <alpha(2)> in N logarithmic steps (set by <steps>). The Tikhonov parameter is specified relative to the largest eigenvalue of M. Consequently, setting alpha higher than 1 strongly restricts the DOF to move away from their initial value. The regularization is stronger for the weak eigenmodes of the problem, which is why it is so well suited for DIC problems. It tends to limit the motion of DOFs which can not be identified anyway. However, it also moves the solution away from the true solution, therefore well defined problems are best solved with low alpha values (e.g. 1e-6). Setting the number of steps to zero causes Tikhonov regularization to be disabled.'...
''...
''...
'Section 7. Correlate'...
'--------------------------------------'...
'<Go> start the correlation process. '...
''...
'<Stop> abort the correlation process, this may take some time before the tool evaluates this action, please be patient. '...
''...
'The <Continue> button sets the current iteration to be converged (regardless of the value of the convergence criteria) and continues to the next CG step.'...
''...
'The status view below the <stop> button shows the correlation progress. In particular it shows which CG step of which increment is currently running, and the value of the convergence parameter for each iteration. Note the "( x/ x) free dof" line, which shows if all DOF are actually used. If due to masking certain basis functions have become ill-supported (less than 20% unmasked) then they are excluded from the DIC algorithm and locked at their initial value.'...
''...
'If the <Live view level> is set to "2", the residual field will be plotted for each iteration. The residual field should reduce to (close to) zero. If a particular area in the ROI has difficulty converging consider adding an initial guess point to that location, or increase the size of the basis functions supported by that area.'...
''...
'Pressing <Go> again starts the correlation process again, however, increments which are complete are not correlated again. To re-correlate an increment it needs to be cleared using the <Clear Inc.> button, which clears the currently visible increment, or the <Clear All> button. Alternatively, use <Restart Inc.> or <Restart All> to the request reusing the previous correlation result as an initial guess to the current correlation, as if it was computed in a previous coarse grain step.'...
''...
'The <to iguess> button takes a complete correlation and converts it to initial guess points, which are added to the existing initial guess points. Typically, it is better to remove the old initial guess points before pressing <to iguess>. After pressing <to iguess> a dialog appears asking to specify how many initial guess points need to be stored (as an x-y grid). It is advised to make sure that the number of initial guess points is more than twice the number of DOF in the coarsest basis. This is to make sure the conversion from initial guess points to DOF is well defined. Otherwise, the tool will approximate the DOF by first fitting a polynomial of adequate order (max 4th order) through the initial guess points.'...
''...
''...
'Section 8. Results'...
'--------------------------------------'...
'The first time this section is selected, the strain fields are computed after selection of the strain definition. Selecting a different strain definition will cause the strain fields to be re-computed.'...
''...
'Three levels of plotting are defined: "Pattern", "Overlay", and "Vectors", which can be plotted simultaneously.'...
''...
'The lowest level is the "Pattern", selecting "f", "gtilde", or "r" (residual) will show the other layers in the undeformed configuration. Selecting "g" will deform the other fields such that they follow the material points. The slider bar below the pattern allows softening of the pattern, i.e. reducing the intensity.'...
''...
'The intermediate level "Overlay", allows plotting of a scalar field, for instance the x-displacement, where the scalar values are displayed as a color. '...
'- U1, and U2 are the x- and y-displacement respectively,'...
'- when using Quasi-3D DIC, U3 is the z-displacement,'...
'- when in 2D, U3, U4, U5 and U6 are constant, linear, quadratic and cubic corrections in'...
'  brightness respectively (with respect to the intensity)'...
'- q is the combined result of all brightness corrections '...
'  q = U3 + U4*gt + U5*gt^2 + U6*gt^3'...
'- r is the residual'...
'The overlay slider bar controls the alpha value of the overlay (0 = transparent, 1 = opaque). On the right side of the colorbar (in the figure panel) are two edit boxes, these allow manual selection of the color limits. If both are 0, then the color limits are automatically scaled.'...
''...
'The top most plot level "Vectors" allows displaying of vector fields, for instance the displacement vectors. The "strain (min,maj)" option shows the major strains and minor strains in their corresponding directions. The slider scales the arrows, if set to 1, then their lengths is to scale.'...
''...
'On the bottom the pixel size (and unit) can be given, such that the displacements are shown in physical quantities. For the strain measures this pixel calibration is not relevant.'...
''...
'The <strain definitions> below are all 2D, switching to a different definition causes the strain tensors to be computed for each pixel in the ROI. Especially the "logarithmic strain" and "Euler-Almansi strain" are slow since they require a tensor spectral decomposition per pixel. The "membrane strain" is for Quasi 3D use, where the 2D strain along the membrane tangent plane is given using a small deformation formulation.'...
''...
'The <Save> button will copy the current figure window and ask where to save the .png file. The copied window is very similar to a normal matlab figure window, thus normal matlab operations are possible. For instance note the "Plot Tools" icon in the top toolbar, which allows changing of almost anything in the figure.'...
''...
'The <Cross-section> button is not implemented yet.'...
''...
'If Arithmetic on these fields is desired the tool can be tricked. In section 9, there are two buttons which allow you to put all the data to the matlab workspace <Guidata to D> and read the (modified) data structure back to the gui <D to Guidata>.'...
''...
'Example:'...
'    <Guidata to D>'...
'    >> D.res(1).Eeq = sqrt(D.res(1).Emin.^2 + D.res(1).Emaj.^2);'...
'    <D to Guidata>'...
''...
'this example reads the data from the gui, computes the equivalent strain field for all pixels for the first increment, and writes the data back to the gui. If the gui now updates its figure (for instance by re-selecting "strain eq." in the Overlay) then the new equivalent strain is shown.'...
''...
''...
'Section 9. Info'...
'--------------------------------------'...
'This Section houses the <Info>, <Status>, <Help> panels. The first to will be filled by using the tool. The data is also stored in D.gui.stat and D.gui.info.'...
''...
'========================================================================'...
' General Tips and Tricks'...
'========================================================================'...
''...
'Tip 1: Initial Guess'...
'--------------------------------------'...
'Spend time on the initial guess, having a good initial guess will make a big difference during correlation.'...
''...
'Tip 2: Coarse basis functions'...
'--------------------------------------'...
'Use full support basis functions (polynomial, harmonic, zernike) for the first two CG steps. These functions are more robust, but have problems capturing the fine details in the deformation fields. The polynomials seem to work well up to 5th order, above, the harmonic basis is recommended.'...
''...
'Tip 3: Use B-Splines'...
'--------------------------------------'...
'Especially the second order B-Splines seem to be a good choice for the final coarse grain step.'...
''...
'Tip 4: Correlate, Learn, Improve, Correlate'...
'--------------------------------------'...
'Start with a quick correlation with simple basis functions and coarse grained images, this makes the correlation fast, and allows identification of the problem areas. Add initial guess points in those problem areas, and try again. The results will not be accurate but will provide an excellent initial guess for the final high detail correlation.'...
''...
'Tip 5: All data is in D'...
'--------------------------------------'...
'Have a look at the guidata structure in the tool <guidata to D>. This structure contains all the data you see in the tool. For instance all the correlation results are stored in "D.cor(inc)" where inc is the increment number. The extra result fields, such as strain are stored in "D.res(inc)", etc.'...
''...
''...
''...
'========================================================================'...
' Back matter'...
'========================================================================'...
''...
'References'...
'--------------------------------------'...
'[1] J. Neggers, J.P.M. Hoefnagels, F. Hild, S.G. Roux and M.G.D. Geers, Direct stress-strain measurements from bulged membranes using topography image correlation, Experimental Mechanics, 2014'...
'[2] Hild F, Roux S. Comparison of Local and Global Approaches to Digital Image Correlation. Experimental Mechanics. 2012 13;52(9):1503–19.'...
'[3] J. Neggers, B. Blaysat, J.P.M. Hoefnagels, M.G.D. Geers, On image gradients in digital image correlation, Int. J. Numer. Meth. Engng, 2015 doi:10.1002/nme.4971'...
''...
''...
'Changelog'...
'--------------------------------------'...
'version 1.00, JN: compiled version of 0.29'...
'version 0.28, JN: the last CG step is now mandatory for more robust correlation, fixed a bug in prev. inc., fixed a bug in updated-lagrange, optimized the plotting code for speed, implemented a new selection tool to replace imrect'...
'version 0.26, JN: improved the dof conversion: masking for pprev, rank test for piguess '...
'version 0.25, JN: added the du convergence parameter'...
'version 0.24, JN: remember save/load locations'...
'version 0.23, JN: added the zernike polynomial basis, correlation reuse and a new defaults implementation'...
'version 0.22, JN: fixed a fundamental flaw in the brightness relaxation (and added up to cubic relaxation)'...
'version 0.21, JN: corrected the green-lagrange strain, and some minor bugs'...
'version 0.20, JN: bug fix release, also updated all the popup messages to better support a multimonitor setup'...
'version 0.19, JN: changed the basis section to allow an unlimited number of basis sets. Added a fourth (extra fine) coarse grain step. Added automatic gradient selection, changed the correlation length algorithm'...
'version 0.18, JN: bugfix release, also reduced the size of save files'...
'version 0.17, JN: improved file management to allow adding and removing of images at any stage in the tool. basicgdic now accepts two optional variables, basicgdic(inputfile,outputfile) to allow it to run headless.'...
'version 0.16, JN: added T6, Q4 and Q8 element types'...
'version 0.14, JN: added autosave feature, modified bestit feature.'...
'version 0.12, JN: bug fixes (many found by Johan)'...
'version 0.10, JN: first beta, all features complete, '...
'version 0.07, JN: alpha state, many changes'...
'version 0.01, JN: initial version, incomplete'...
''...
'Contact Information'...
'--------------------------------------'...
'janneggers@gmail.com'...
''...
''...
'Information about the Code Structure'...
'--------------------------------------'...
''...
'The most interesting files are:'...
'basicgdic.m --> contains most of the functions directly connected to a control'...
'lib_bgdic/basicgdic_uicontrols.m --> creates all the controls (buttons,axes,etc.)'...
'lib_bgdic/corsequence.m --> Controls the correlation sequence (pre-processing)'...
'lib_bgdic/correlate.m --> Correlates two images'...
'lib_bgdic/rescompute.m --> Compute the results (strain fields)'...
''...
'Related to the GUI:'...
'The main window of the graphical user interface (GUI) is a matlab figure window. The figure window can be identified programmatically by it''s handle. If it is the first window you open after starting matlab, the handle will be 1. Throughout the code this handle is stored in the variable H. Any matlab figure window can also contain data. This data can be stored and loaded using the guidata function. This mechanism is used throughout the code to transfer the data from one function to the next. Therefore, most functions start with, <D = guidata(H)> which reads the data from figure window, and end with, <guidata(H,D)> to store the data back to the figure window. Functions which do not have the <D = guidata(H)> line simply do not require any data, and functions which do not have the <guidata(H,D)> line simply do not change any data.'...
''...
'Each control (button/axes/etc.) has a list of properties connected to it, two of which are interesting: ''Tag'' and ''Call''. The Tag of a control is used to generate an automatic structure of handles using the <S = guihandles(H)> command. These handles can then be used to modify all properties of the control such as the ''Value'' or the ''String''. The Call property specifies a function which is called when that particular control is activated (e.g. pressing a button). All the (smaller) callback functions are inside the basicgdic.m file, and are grouped per section, with the general functions listed before the first section.'...
''...
'The m-file lib_bgdic/basicgdic_uicontrols.m is called in the beginning of the basicgidc.m function. In this m-file all the panels, buttons, axes, etc. are defined. They are organized per section, with the general definitions in the top (for the panels and the sections itself). So if you wish to find the properties of a certain button, look for the big comment indication it''s section, and then the particular button is nearby. This GUI contains three panels, which seem to switch content. This is not actually true, the way this is programmed is that concurrently 9 versions of the control panel and figure panel exist, but 8 of them are invisible. Pressing the buttons in the section selection panel only changes the visibility state such that the correct panel is now visible. Every object which is a child (or grandchild) of a panel is automatically set to the same visibility state as the panel. Therefore, setting the panel to invisible causes every object inside the panel to become invisible if it''s parent property is set correctly. Consequently, if you wish to add a button to a panel you need to set the parent property correctly such that it''s visibility state is handled correctly.'...
''...
'Have a look at the lib_btgdic/plot<something>.m functions to see how the figures in the various figure panels are created.'...
''...
'Related to DIC:'...
'A DIC code, which correlates the displacement field between two images is really very simple and rather short. Including the code to deal with all the options in this tool, the lib_bgdic/correlate.m function is still less than 1000 lines. Strip it to the basics and you can write a GDIC code in less than 20 lines. However, more work needs to be done before a correlation can start. How to transfer the results from one correlation to the next, considering that the previous correlation can be a different coarse grain step, or a different increment, or the first increment, or etc? All this work is handled by the lib_bgdic/corsequence.m file, with some helper functions. It processes the previous results and the settings for the next correlation and prepares the new basis functions and the initial values of the corresponding DOF. Afterwards, some post-processing is required (computing the strain fields) which is performed in the lib_bgdic/rescompute.m file.'...
''...
'Please do not hesitate to email me if you plan to modify this program and need some help.'...
''...
''...
'Copyright'...
'--------------------------------------'...
''...
''...
'Copyright 2015 Jan Neggers'...
''...
'This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.'...
''...
'This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.'...
''...
''...
''...
''...
''...
''...
};

