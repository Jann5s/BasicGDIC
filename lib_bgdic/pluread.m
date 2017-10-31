function plu = pluread(varargin)
% plu = pluread(filename)
%
% This function reads a .plu file from the Sensofar optical profilomter it
% returns a structure with information about the measurement and the data
% itself (plu.Z), lost pixels are set to NaN. The data is stored in the
% .plu file as (32bit) single precision floating point, however, pluread
% converts this to double precision with the dimension of micrometers.
%
% The plu software allows for storing 4 types of data (Topography, Image,
% Profile or Single point, and, Multiple Profile). Currently only the first
% two are supported (implementing the others should be easy).
%
% Example:
% plu = pluread('somefile.plu')
% x = plu.x ;
% y = plu.y ;
% Z = plu.Z ;
% surf(x,y,Z,'EdgeAlpha',0.1,'FaceAlpha',0.5)
%
% The figure created by "surf" is nice but very slow, a quicker option
% would be:
% imagesc(x,y,Z)


% Changelog:
% version 1.2 19-oct-12
% Lambert Bergers; added code to extract image from "conventional image"
%                  type of plu-files.
% version 1.1 5-Jun-12
% Jan Neggers; created pcode file
%
% version 1.0 13-Dec-10
% Jan Neggers; cleaned the output
%
% version 0.94 02-jun-09
% Jan Neggers; added plu.x and plu.y to the structure from the old example
%              plu.x = plu.mppx*(0:plu.N(1)-1) ;
%              plu.y = plu.mppy*(0:plu.M(1)-1) ;
%              and fixed the transpose booboo in the old example
%
% version 0.93 09-apr-09
% Norman Delhey; fixed a bug when no stack images are present in the .plu
%
% version 0.92 04-mar-09
% Jan Neggers; added "filedate" and "filedatenum", which are based on the
%              modification date of the file.
%
% version 0.91 16-feb-09
% Jan Neggers; fixed some minor bugs
%
% version 0.9  08-feb-09
% Jan Neggers; first release

if nargin == 0
    help pluread
    fprintf('\n\n\n')
    return
else
    filename = varargin{1};
end
if strcmpi(filename,'help')
    help pluread
    return
elseif strcmpi(filename,'changelog')
    return
end


try
    
    % fix the extention, 
    if isempty(regexpi(filename,'\.plu$','match'))
        filename = [filename '.plu'];
    end
    
    % storing the filename
    plu.file = filename;
    
    % preallocating fields in the structure (to define the order in which they show)
    plu.filename    = '';
    plu.filesize    = 0 ;
    plu.filedate    = '';
    plu.filedatenum = '';
    
    % Filling the tables with constants
    plu.CTables.C1.Max_Long_Data = 128;
    plu.CTables.C1.LostPixels    = 1000001;
    
    % C2 Measurement Types
    plu.CTables.C2{1}  = 'Confocal Image';
    plu.CTables.C2{2}  = 'Profile';
    plu.CTables.C2{3}  = 'Multiple profile';
    plu.CTables.C2{4}  = 'Topography';
    plu.CTables.C2{5}  = 'Coordinate';
    plu.CTables.C2{6}  = 'Single Point Thickness';
    plu.CTables.C2{7}  = 'Custom Application';
    
    % C3 Acquisition Algorithms
    plu.CTables.C3{1}  = 'Confocal Intensity';
    plu.CTables.C3{2}  = 'Confocal Gradient';
    plu.CTables.C3{3}  = 'Interferometric PSI';
    plu.CTables.C3{4}  = 'Interferometric VSI';
    plu.CTables.C3{5}  = 'Interferometric ePSI';
    plu.CTables.C3{6}  = 'Confocal thickness';
    plu.CTables.C3{7}  = 'Interferometric thickness';
    
    % C4 Acquisition Methods
    plu.CTables.C4{1}  = 'Conventional Image';
    plu.CTables.C4{2}  = 'Confocal Image';
    plu.CTables.C4{3}  = 'Single Profle';
    plu.CTables.C4{4}  = 'Extended Profle';
    plu.CTables.C4{5}  = 'Topography';
    plu.CTables.C4{6}  = 'Extended Topography';
    plu.CTables.C4{7}  = 'Multiple Profile';
    plu.CTables.C4{8}  = 'Extended Multiple Profile';
    plu.CTables.C4{9}  = 'Coordinate';
    plu.CTables.C4{10} = 'Custom';
    
    %C5 Objectives
    plu.CTables.C5{1}  = 'Nikon SLWD 10x';
    plu.CTables.C5{2}  = 'Nikon SLWD 20x';
    plu.CTables.C5{3}  = 'Nikon SLWD 50x';
    plu.CTables.C5{4}  = 'Nikon SLWD 100x';
    plu.CTables.C5{5}  = 'Nikon EPI 20x';
    plu.CTables.C5{6}  = 'Nikon EPI 50x';
    plu.CTables.C5{7}  = 'Nikon EPI 10x';
    plu.CTables.C5{8}  = 'Nikon EPI 100x';
    plu.CTables.C5{9}  = 'Nikon ELWD 10x';
    plu.CTables.C5{10} = 'Nikon ELWD 20x';
    plu.CTables.C5{11} = 'Nikon ELWD 50x';
    plu.CTables.C5{12} = 'Nikon ELWD 100x';
    plu.CTables.C5{13} = 'Nikon 2.5x TI';
    plu.CTables.C5{14} = 'Nikon 5x TI';
    plu.CTables.C5{15} = 'Nikon 10x DI';
    plu.CTables.C5{16} = 'Nikon 20x DI';
    plu.CTables.C5{17} = 'Nikon 40x DI';
    plu.CTables.C5{18} = 'Nikon EPI 5x';
    plu.CTables.C5{19} = 'Nikon EPI 150x';
    plu.CTables.C5{20} = 'Unknown';
    
    % C6 Field of view areas
    plu.CTables.C6{1}  = '128 x 128';
    plu.CTables.C6{2}  = '256 x 256';
    plu.CTables.C6{3}  = '512 x 512';
    plu.CTables.C6{4}  = 'CCDrows x CCDcols';
    plu.CTables.C6{5}  = '256 x CCDcols';
    plu.CTables.C6{6}  = '128 x CCDcols';
    
    % C7 File format version
    plu.CTables.C7{1}  = '2000';
    plu.CTables.C7{2}  = '2006';
    
    % C8 Hardware configuration
    plu.CTables.C8{1}  = 'PLU';
    plu.CTables.C8{2}  = 'PLU 2300 XGA (2003)';
    plu.CTables.C8{3}  = 'PLU 2300 XGA T5 (2004)';
    plu.CTables.C8{4}  = 'PLU 2300 SXGA (2005)';
    plu.CTables.C8{5}  = 'PLU 3300';
    
    % Reading the file details
    D = dir(plu.file);
    if length(D) ~= 1
        error('Specify a file in stead of a directory as input.');
    end
    plu.filesize    = D.bytes;
    plu.filename    = D.name;
    plu.filedate    = D.date;
    plu.filedatenum = D.datenum;
    
    % Opening the file
    [fid fopenmessage] = fopen(plu.file);
    if fid == -1
        error(fopenmessage)
    end
    
    % Reading the machineformat from the file
    [filename,permission,machineformat,encoding] = fopen(fid);
    
    % Setting some counters
    plu.bytes.A1 = 0; plu.bytes.A2 = 0; plu.bytes.A3 = 0;
    
    %%% Reading Area 1 - The Header
    %===============================
    
    % The date
    count     = 128; % tell fread to read 128 blocks
    precision = 'schar=>char'; % tell fread what type of blocks
    
    plu.date = fread(fid, count, precision, 0, machineformat).';
    plu.date = deblank(plu.date); % cleaning the string of blanks
    plu.bytes.A1 = plu.bytes.A1 + count*1; % counting the processed bytes
    
    % The date in numerical format (can't figure out which format)
    count     = 1;
    precision = 'int32=>double'; %???? not working (check datevec(plu.datenum))
    plu.datenum = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    % comment
    count     = 256;
    precision = 'schar=>char';
    plu.CTables.Comment = fread(fid, count, precision, 0, machineformat).';
    plu.CTables.Comment = deblank(plu.CTables.Comment);
    plu.bytes.A1 = plu.bytes.A1 + count*1;
    
    % XY callibration (10 blocks in total)
    count     = 3;
    precision = 'uint32=>double';
    temp = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    plu.N         = temp(1);  % Y axis # of data points
    plu.M         = temp(2);  % X axis # of data points
    plu.N_CutAxis = temp(3);  % Cut axis # of points
    
    count     = 7;
    precision = 'float32=>double';
    temp = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    plu.dy_multip    = temp(1);  % Y axis seperation between profiles of a multiple profile (um)
    plu.mppx         = temp(2);  % X axis um per pixel
    plu.mppy         = temp(3);  % Y axis um per pixel
    plu.mpp_CutAxis  = temp(6);  % Cut axis um per pixel (reordered before the x0)
    plu.x0           = temp(4);  % X axis coordinate of first data point
    plu.y0           = temp(5);  % Y axis coordinate of first data point
    plu.z0           = temp(7);  % z axis coordinate of mean value
    
    % Measurement configuration
    count     = 10;
    precision = 'uint32=>double';
    temp = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    % temp 1 2 4 5 select directly from the tables above except that their
    % numbering starts at 0 and the cell index starts at 1 (hence the +1)
    plu.MeasurementType         = plu.CTables.C2{temp(1)+1};
    plu.AcquisitionAlgorithm    = plu.CTables.C3{temp(2)+1};
    % temp 3 gives only 0 or 1, requiring a bit more complex selecting based on
    % temp 2 also.
    plu.AcquisitionMethod       = [ 0 2 6 4 8 9 9 ];
    plu.AcquisitionMethod       = plu.CTables.C4{plu.AcquisitionMethod(temp(2)+1)+temp(3)};
    plu.Objective               = plu.CTables.C5{temp(4)+1};
    plu.FieldOfView             = plu.CTables.C6{temp(5)+1};
    plu.N_FOV                   = temp(6);
    plu.M_FOV                   = temp(7);
    plu.N(2)                    = temp(8); % info already obtained from XYcal ??
    plu.M(2)                    = temp(9); % info already obtained from XYcal ??
    plu.NumberOfFOVs            = temp(10);
    
    count     = 1;
    precision = 'double=>double';
    plu.z_incr = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*8;
    
    count     = 1;
    precision = 'float=>double';
    plu.z_range = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    count     = 1;
    precision = 'uint32=>double';
    plu.z_planes = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    count     = 1;
    precision = 'uint32=>double';
    plu.AcquisitionThreshold = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    count     = 1;
    precision = 'int8=>double';
    plu.Restore = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*1;
    
    count     = 5;
    precision = 'uchar=>uint8';
    temp = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*1;
    
    plu.NumberOfLayers        = temp(1);
    plu.FileFormatVersion     = temp(2);
    plu.HardwareConfiguration = plu.CTables.C8{temp(3)+1};
    plu.NumberOfStackImages   = temp(4);
    plu.NotUsed01             = temp(5);
    
    count     = 1;
    precision = 'int32=>double';
    plu.DecimationFactor = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*4;
    
    % adding vectors for the x and y axis
    plu.x = plu.mppx*(0:plu.M(1)-1) ;
    plu.y = plu.mppy*(0:plu.N(1)-1) ;
    
    % skipping 2 (unused?) bytes
    count     = 2;
    precision = 'uchar=>double';
    temp = fread(fid, count, precision, 0, machineformat);
    plu.bytes.A1 = plu.bytes.A1 + count*1;
    
    
    %%% Reading Area 2 - The Data
    %===============================
    % There are 4 (+1 misc) types of plu files, the raw data configuration in the file
    % depends on the type listed in plu.MeasurementType
    
    if strcmp(plu.MeasurementType,plu.CTables.C2{4})
        % Configuration 1: Topography (tested)
        
        % reading the number of rows and columns
        count     = 2;
        precision = 'uint32=>double';
        plu.RawSize = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
        % reading the data into one matrix
        count     = [plu.RawSize(2) plu.RawSize(1)];
        precision = 'float32=>double';
        plu.Z = fread(fid, count, precision, 0, machineformat).';
        plu.Z(plu.Z==plu.CTables.C1.LostPixels) = NaN;
        plu.bytes.A2 = plu.bytes.A2 + prod(count)*4;
        
        % Max value of the data
        count     = 1;
        precision = 'float32=>double';
        plu.Max = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
        % Min value of the data
        count     = 1;
        precision = 'float32=>double';
        plu.Min = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
        % Reading extra added images
        if plu.NumberOfStackImages ~= 0
            for i = 1:plu.NumberOfStackImages
                clear A R G B
                
                count = prod(plu.RawSize)*3;
                precision = 'uint8=>uint8';
                
                % reading the image into one big column
                StackImage = fread(fid, count, precision, 0, machineformat).';
                plu.bytes.A2 = plu.bytes.A2 + prod(count)*1;
                
                if length(StackImage) == 3*plu.RawSize(2)*plu.RawSize(1)
                    % reshaping it to a (3,N,M) matrix
                    A = reshape(StackImage,3,plu.RawSize(2),plu.RawSize(1));
                    % isolating the R,G,B matrices
                    R = double(squeeze(A(1,:,:))).' ./ 255;
                    G = double(squeeze(A(2,:,:))).' ./ 255;
                    B = double(squeeze(A(2,:,:))).' ./ 255;
                    
                    % converting to 16bit grayscale
                    plu.Image{i} = uint16((0.299*R + 0.587*G + 0.114*B)*65535) ;
                end
                
            end
        end
    elseif strcmp(plu.MeasurementType,plu.CTables.C2{1})
        % Configuration 2: Image
        
        % reading the number of rows and columns
        count     = 2;
        precision = 'uint32=>double';
        plu.RawSize = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
        % reading the data into one matrix
        count     = [plu.RawSize(2) plu.RawSize(1)];
        precision = 'float32=>double';
        plu.Z = fread(fid, count, precision, 0, machineformat).';
        plu.Z(plu.Z==plu.CTables.C1.LostPixels) = NaN;
        plu.bytes.A2 = plu.bytes.A2 + prod(count)*4;
        % Note, Sensofar actually stores uint8 data as single precision
        % float
        
        % Max value of the data
        count     = 1;
        precision = 'float32=>double';
        plu.Max = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
        % Min value of the data
        count     = 1;
        precision = 'float32=>double';
        plu.Min = fread(fid, count, precision, 0, machineformat).';
        plu.bytes.A2 = plu.bytes.A2 + count*4;
        
    elseif strcmp(plu.MeasurementType,plu.CTables.C2{2}) || strcmp(plu.MeasurementType,plu.CTables.C2{6})
        % Configuration 3: Profile or Single point thickness
        error('This part of the pluread is not ready yet, DataType Profile or Single point')
    elseif strcmp(plu.MeasurementType,plu.CTables.C2{3})
        % Configuration 4: Multiple Profile
        error('This part of the pluread is not ready yet, DataType Multiple Profile')
    else
        % Configuration 5: Other measurement type
        error('Unknown Data type in this plu file, try the plu software')
    end
    
    % Adding all processed bytes for debugging/reviewing
    plu.bytes.processed = plu.bytes.A1 + plu.bytes.A2 + plu.bytes.A3;
    
    % closing the file
    fclose(fid);
    
catch ME
    % If an error occured, show the error
    disp(ME)
    error(ME.message)
    if exist('fid','var')
        % and if the file is still open close it
        fclose(fid);
    end
end

% ======================================
% Cleanup the output
% ======================================

% imporatant output
P.file      = plu.file;
P.date      = plu.date;
P.pixelsize = [ plu.mppx plu.mppy ];
P.x         = plu.x;
P.y         = plu.y;
P.Z         = plu.Z;
P.x0        = plu.x0;
P.y0        = plu.y0;
P.z0        = plu.z0;

% File info
F.filename    = plu.filename;
F.filesize    = plu.filesize;
F.filedate    = plu.filedate;
F.filedatenum = plu.filedatenum;

% Measurement info
M.N_CutAxis             = plu.N_CutAxis;
M.mpp_CutAxis           = plu.mpp_CutAxis;
M.dy_multip             = plu.dy_multip;
M.MeasurementType       = plu.MeasurementType;
M.AcquisitionAlgorithm  = plu.AcquisitionAlgorithm;
M.AcquisitionMethod     = plu.AcquisitionMethod;
M.Objective             = plu.Objective;
M.FieldOfView           = plu.FieldOfView;
M.FOV                   = [plu.N_FOV plu.M_FOV];
M.NumberOfFOVs          = plu.NumberOfFOVs;
M.z_incr                = plu.z_incr;
M.z_range               = plu.z_range;
M.z_planes              = plu.z_planes;
M.AcquisitionThreshold  = plu.AcquisitionThreshold;
M.Restore               = plu.Restore;
M.HardwareConfiguration = plu.HardwareConfiguration;
M.NumberOfStackImages   = plu.NumberOfStackImages;
M.DecimationFactor      = plu.DecimationFactor;

% Extra info
E.N           = plu.N;
E.M           = plu.M;
E.bytes       = plu.bytes;
E.CTables     = plu.CTables;
E.datenum     = plu.datenum;
E.RawSize     = plu.RawSize;
if isfield(plu,'Image')
    E.Image       = plu.Image;
end

% restore the output
clear plu

plu = P;
plu.FileInfo = F;
plu.MeasurementInfo = M;
plu.ExtraInfo = E;