function scal = DisplayUndistortionCSV(calibinfilename, caliboutfilename, screenid)
% scal = DisplayUndistortionCSV(calibinfilename [, caliboutfilename][, screenid])
%
% Geometric display calibration procedure for undistortion of distorted
% displays. Needs graphics hardware with basic support for the PTB imaging
% pipeline (see below).
%
% Many display devices, e.g., video beamers and most CRT displays cause
% some amount of spatial distortion to your visual stimuli during display.
% Video projection onto curved screens will cause an especially large
% amount of distortion.
%
% Psychtoolbox can "undistort" your visual stimuli for you: At stimulus
% onset time, PTB applies a geometric warping transformation to your
% stimulus which is meant to counteract or cancel out the geometric
% distortion caused by your display device. If both, PTB's warp transform
% and the implicit distortion transform of the display match, your stimulus
% will show up undistorted on the display device.
%
% For this to work, PTB needs two things:
%
% 1. Recent graphics hardware with support for the PTB imaging pipeline:
% See our Wiki for recommendations. However, all ATI cards starting with
% Radeon 9500 and all NVidia cards of type GeForce-FX5200 and later, as
% well as the Intel-GMA 950 and later should be able to do it, although
% more recent cards will have a higher performance.
%
% 2. A calibration file that defines the warp transformation to apply. Your
% experiment script will load that file into Screen's "warp engine" at the
% beginning of your experiment and Screen() will automatically apply that
% warping to each stimulus image before output.
%
% DisplayUndistortionCSV defines a continous mapping (x', y') = f(x, y)
% from uncorrected input pixel locations (x,y) in your stimulus image to
% output locations (x', y') on your display. This mapping is defined by a
% linear mesh of quadrilaterals, as read from the given input file. The
% file format is compatible with the Warp-API of NVidia's proprietary
% graphics drivers for some high-end Quadro GPUs on Windows-7 and later.
%
%
% How to use:
% -----------
%
%
% Execute the function with the following parameters:
%
%
% `calibinfilename` is the path and filename of an existing calibration
% file. This must be a ASCII file with a one-line header, describing the
% meaning of the data columns, followed by at least 4 lines with numeric
% calibration data. Each line defines a data row, which has six columns,
% delimited with a semicolon ';'. The columns of a row have the following
% meaning: vx;vy;tx;ty;c;r
%
% vx = Output mesh vertex x position.
% vy = Output mesh vertex y position.
% tx = Input framebuffer corresponding x position.
% ty = Input framebuffer corresponding y position.
% c  = Column index of the vertex.
% r  = Row index of the vertex.
%
% Each row describes one vertex of a mesh of connected quadrilaterals
% (GL_QUADS). The mesh is basically a rectilinear grid, whose "grid line"
% intersections are the locations of the vertices. The mesh gets distorted
% by moving around those intersections, thereby forming the actual warp
% mesh.
%
% Pixel color values at location (tx,ty) of the input framebuffer are
% projected to location (vx,vy) of the warped output framebuffer. Quads are
% defined by 4 vertices which define their corners. Pixel sampling and
% output locations inside quads are bilinearly interpolated for a smooth
% warp grid. The origin of the calibration coordinate system is the
% top-left framebuffer / screen corner (0,0), x-axis points to the right,
% y-axis points downwards. All coordinates are normalized to a unit square
% screen. The range 0.0 - 1.0 maps to the screen width or height (0.0 =
% Left or Top, 1.0 = Right or bottom) and is subject to rescaling for a
% specific displays real resolution. Values outside the 0-1 range are
% allowed. The corresponding mesh locations will be clipped to screen
% corners during warping.
%
% The file format described here is the format expected by NVidia's
% Warp-API, as implemented in some binary NVidia graphics drivers for some
% high-end Quadro graphics cards on Windows-7 and later. This allows to
% reuse the file for both NVidia proprietary distortion correction and for
% Psychtoolbox.
%
%
% `caliboutfilename` Name of the file to which calibration results should
% be stored. If no name is provided, the file will be stored inside the
% 'GeometryCalibration' subfolder of your psychtoolbox configuration
% directory (path is PsychToolboxConfigDir('GeometryCalibration'). The
% filename will contain the screenid of the display that was calibrated.
%
%
% `screenid` screen id of the target display for calibration. The parameter
% is optional, defaults to zero, and is only used to generate the default
% filename for the output file.
%
%
% This script will print out a little snippet of code that you can paste
% and include into your experiment script - That will automatically load
% the calibration result file and apply the proper undistortion operation.
%

% History:
%
% 07/27/12  Written (MK).

% Child protection:
AssertOpenGL;

if ~exist('screenid', 'var') || isempty(screenid)
    screenid = 0;
end

if ~exist('calibinfilename', 'var') || isempty(calibinfilename) || ~ischar(calibinfilename)
    error('You must provide the namestring of a valid calibration input file as parameter "calibinfilename"!');
end

% Check for such a file:
if ~exist(calibinfilename, 'file')
    error('Calibration input file %s does not exist or is inaccessible!', calibinfilename);
end

% Ok, valid file provided. Read it in:
try
    % Parse ASCII text file. First line is header - skip it, start at row
    % index 1 instead. Start with first column (0). Delimiter is a
    % semicolon. Matrix m contains all values:
    m = dlmread(calibinfilename, ';', 1, 0);
catch %#ok<*CTCH>
    error('Failed to load input calibration file %s for some reason.', calibinfilename);
end

fprintf('\nLoaded calibration data from ASCII file %s.\n', calibinfilename);

% At this point, screenid contains the final screenid for the screen to
% calibrate. Assign it to scal struct. This will create 'scal' if it
% doesn't exist yet, or override its screenid in some cases:
scal.screenNumber = screenid;

% Define type of mapping for this calibration method: This is used in the
% CreateDisplayWarp() routine when parsing the calibration file to detect
% the type of undistortion method to use, ie. how to interpret the data in
% the calibration file to setup the calibration.
warptype = 'CSVDisplayList'; %#ok<NASGU>

n = size(m, 1);
if n < 4 || size(m, 2) ~= 6
    error('ASCII input file must have at least 4 rows with six columns of data in each row, to define at least one quadrilateral.')
end

% Find number of vertex rows and columns for calibration mesh:
cols = max(m(:, 5)) + 1; % Column indices are zero-based, so add 1 for real count.
rows = max(m(:, 6)) + 1; % Row indices are zero-based, so add 1 for real count.

if n ~= (cols * rows)
    error('Mismatch between number of data lines in file and maximum row and column indices in ASCII file! File corrupt?!?');
end

% Build 2D matrices: rows x cols per plane, 2 planes for x and y
% components:

% Vertex coordinates of the rendered output meshes quad vertices:
scal.vcoords = zeros(rows, cols, 2);

% Corresponding texture coordinates for sourcing from input framebuffer:
scal.tcoords = zeros(rows, cols, 2);

% Parse the 'm' matrix and fill scal.vcoords and scal.tcoords:
for i = 1:n
    tc = m(i, 5) + 1; % Target column in matrices.
    tr = m(i, 6) + 1; % Target row in matrices.

    scal.vcoords(tr, tc, 1) = m(i, 1); % vertex x.
    scal.vcoords(tr, tc, 2) = m(i, 2); % vertex y.

    scal.tcoords(tr, tc, 1) = m(i, 3); % texture x [u].
    scal.tcoords(tr, tc, 2) = m(i, 4); % vertex y [v].
end

% 'scal' contains the final results of calibration. Write it out to
% calibfile for later use by the runtime routines:

% Check if name for calibration result file is provided:
if ~exist('caliboutfilename', 'var')
    caliboutfilename = [];
end

if isempty(caliboutfilename)
    % Nope: Assign default name - Store in dedicated subfolder of users PTB
    % config dir, with a well defined name that also encodes the screenid
    % and resolution for which to calibrate:
    caliboutfilename = [ PsychtoolboxConfigDir('GeometryCalibration') 'CSVCalibdata' sprintf('_%i', screenid) '.mat'];
    fprintf('\nNo name for calibration file provided. Using default name and location...\n');
end

% Print name of calibfile and check for existence of file:
fprintf('Name of calibration result file: %s\n\n', caliboutfilename);
if exist(caliboutfilename, 'file')
    answer = input('This file already exists. Overwrite it [y/n]? ','s');
    if ~strcmpi(answer, 'y')
        fprintf('\n\nCalibration aborted. Please choose a different name for calibration result file.\n\n');
        return;
    end
end

% Save all relevant calibration variables to file 'caliboutfilename'. This
% method should work on both, Matlab 6.x, 7.x, ... and GNU/Octave - create
% files that are readable by all runtime environments:
save(caliboutfilename, 'warptype', 'scal', '-mat', '-V6');

fprintf('Import of Calibration file finished :-)\n\n');
fprintf('You can apply the calibration in your experiment script by replacing your \n')
fprintf('win = Screen(''OpenWindow'', ...); command by the following sequence of \n');
fprintf('commands:\n\n');
fprintf('PsychImaging(''PrepareConfiguration'');\n');
fprintf('PsychImaging(''AddTask'', ''LeftView'', ''GeometryCorrection'', ''%s'');\n', caliboutfilename);
fprintf('win = PsychImaging(''OpenWindow'', ...);\n\n');
fprintf('This would apply the calibration to the left-eye display of a stereo setup.\n');
fprintf('Additional options would be ''RightView'' for the right-eye display of a stereo setup,\n');
fprintf('or ''AllViews'' for both views of a stereo setup or the single display of a mono\n');
fprintf('setup.\n\n');
fprintf('The ''GeometryCorrection'' call has a ''debug'' flag as an additional optional parameter.\n');
fprintf('Set it to a non-zero value for diagnostic output at runtime.\n');
fprintf('E.g., PsychImaging(''AddTask'', ''LeftView'', ''GeometryCorrection'', ''%s'', 1);\n', caliboutfilename);
fprintf('would provide some debug output when actually using the calibration at runtime.\n\n\n');

% Done.
return;
