function [scal] = DisplayUndistortionBVL(caliboutfilename, screenid, xnum, ynum, referenceImage, stereomode)
% [scal] = DisplayUndistortionBVL([caliboutfilename][, screenid][, xnum=37][, ynum=27][, referenceImage=None][, stereomode=0])
% [scal] = DisplayUndistortionBVL([caliboutfilename][, calibinfilename])
%
% Geometric display calibration procedure for undistortion of distorted
% displays. Needs graphics hardware with basic support for the PTB imaging
% pipeline (see below).
%
% This code was contributed by the members of the Banks Vision Lab at
% University of California, Berkeley. It is a subset of their internal
% "BVL" library and used since many years for their display devices,
% haploscopes etc., so it should be reasonably bug-free and mature.
%
% Many display devices, e.g., video beamers and most CRT displays cause
% some amount of spatial distortion to your visual stimuli during display.
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
% beginning of your experiment.
%
% This routine allows you to create such a calibration file in an
% interactive procedure: You can change and tweak the transformation until
% it fits your needs, ie. until it nicely undistorts your display. Then the
% corresponding calibration file is saved for later use with that display.
%
% DisplayUndistortionBVL defines a continous mapping (x', y') = f(x, y)
% from uncorrected input pixel locations (x,y) in your stimulus image to
% output locations (x', y') on your display. This mapping is defined by a
% 3rd order, 2-dimensional polynomial that is fitted to the calibration data.
%
%
% How to use:
% 
% 1. Start the script, providing all parameters that you don't want to have
% at default settings (all parameters have reasonable defaults):
%
% 'caliboutfilename' Name of the file to which calibration results should
% be stored. If no name is provided, the file will be stored inside the
% 'GeometryCalibration' subfolder of your psychtoolbox configuration
% directory (path is PsychToolboxConfigDir('GeometryCalibration'). The
% filename will contain the screenid and resolution of the display that was
% calibrated.
%
% 'calibinfilename' is the optional path and filename of an existing
% calibration file that you want to use as a template or starting point for
% your calibration, e.g., from previous calibrations in order to save you
% some setup work. If no 'calibinfilename' is provided or if 'screenid' and
% other parameters are provided instead, the script will start with a
% default rectilinear calibration grid of equally spaced points.
%
% These parameters only have meaning if no 'calibinfilename' was provided:
%
% 'screenid' screen id of the target display for calibration. If
% 'calibinfilename' was provided, the screenid encoded in that file will be
% used. Otherwise, the provided screenid will be used. If this parameter is
% omitted, PTB will use the single screen on a single display setup. On a
% multi-display setup, PTB will ask for the screenid.
%
% 'xnum' and 'ynum' Optional number of horizontal and vertical calibration
% points to use. Their number doesn't affect runtime behaviour of your
% stimulus script, only the quality of your calibration and the amount of
% time you'll need to calibrate. The default settings are xnum = 37 and
% ynum = 27. 'xnum' and 'ynum' should be odd numbers. If you provide an
% even number it will be rounded up to the next odd number.
%
% 'referenceImage' Optional name of an image file in a format supported by
% Matlab/Octaves imread() command: If provided, the image will be loaded
% from filesystem and drawn as a backdrop to the calibration grid. You can
% use this if you want to use this routine not to undistort a physical
% display, but want to undistort an existing image, e.g., create a proper
% calibration file for the ImageUndistortionDemo routine.
%
% 'stereomode' Optional stereo mode for display of calibration. Defaults to
% zero, i.e., monoscopic display.
%
% 2. After startup, the script will display a grid of mostly evenly spaced
% points onscreen. The points will not be perfectly aligned to a grid due
% to the distortion caused by your display. Your job is to tweak and shift
% those points so that they line up to a rectilinear grid as good as
% possible on your display from the viewpoint of your subject.
%
% The best way to do this is to build a real mechanical rectilinear grid of
% fine wires and mount it to your display as a "real world" reference for a
% perfect grid. Then you move the displayed calibration dots so that they
% perfectly align with the crossings of the horizontal- and vertical wires
% of your mechanical reference grid.
%
% Mouse operation:
%
% In 'global mode' (when you see the hair-cross in the center of the
% screen), the mouse buttons do the following:
%
% 'Left mouse' or 'l' key:   Switch to local mode.
% 'Middle mouse' or 'm' key: Change stepsize for parameter adjustments.
% 'Right mouse' or 'r' key:  Change type of parameter to adjust.
%
% Following global parameters can be adjusted in global mode:
% 'Translation' of whole dot field.
% 'Scale' of dot field.
% 'Shear' of dot field.
% 'Trapezoid' correction.
% 'Barrel and pincusion' distortion.
%
% In 'local mode': The location of a selected dot or a selection of dots
% can be changed:
%
% 'Left mouse' + mouse drag: Select a single dot close to mouse cursor, or
%                            draw a bounding rectangle around a region of
%                            dots to select.
%
% 'Middle mouse' or 'm' key: Change stepsize for translation of dot(s).
% 'Right mouse' or 'r' key:  Unselect dot(s) if dot(s) is/are selected.
%                            Switch to 'global mode' if no dot(s) selected.
%
%
% 
% Keys and their meaning:
%
% 'l' 'm' and 'r' buttons are synonyms for left- middle- and right mouse
% buttons. However, its much more convenient to use a three button mouse.
%
% Cursor arrow keys move selected dot(s) or change selected global
% parameters.
%
% 'space' key toggles the display of the online help text.
%
% You finish the calibration and write it into a calibration file by
% pressing the ESCape key. This will end the calibration script. The
% structure 'scal' which contains all calibration results will also be
% returned by this function as optional return argument.
%
% This script will print out a little snippet of code that you can paste
% and include into your experiment script - That will automatically load
% the calibration result file and apply the proper undistortion operation.
%
%

% History:
%
% 02/17/08  Initial (incomplete) proposal, derived from David Hoffmanns
%           example script script4David.m. (MK)
%           This is just an initial working draft: Its basic naming
%           conventions and syntax are similar to the
%           DisplayUndistortionBezier.m calibration script.
%
% 05/01/08  Refinement: Integrated all BVL helper routines into distro, as
%           hidden (private) subfolder of PsychGLImageProcessing, so
%           they're only accessible from this master function.
%
%           Fixed and modified a few of the helpers, moved some init and
%           error checking code etc. to this master function for a
%           "smoother" user experience. Now user can either pass in the
%           name of a previous calibration file (as starting point for a
%           calibration), or create a calibration from scratch. If creating
%           from scratch, the 'screenid' can be specified on multi-display
%           systems. On single display setups, the single screen is used
%           and on multi-display setups without 'screenid', the user is
%           asked.
%
%           The default number of calibration dots is 37 by 27, but user
%           can override via the 'xnum' and 'ynum' arguments.
%
%           Output files are now stored in a subfolder of
%           PsychtoolboxConfigDir, with a name that encodes screenid and
%           resolution by default.
%
% 05/03/13  Add optional stereomode parameter to select other display modes
%           than monoscopic display for calibration. (MK)

% Running on PTB-3? Abort otherwise:
PsychDefaultSetup(1);

% Need to temporarily add the '/private' subfolder to Octave's path, as
% Octave as of V3.0.5 doesn't know about private subfolders yet:
if IsOctave
    warning('As of Octave version 3.2.0, this function does not always produce correct results, due to insufficient support for the griddata() v3 interpolation method!');
    addpath([ fileparts(mfilename('fullpath')) '/private' ]);
end

% Enable unified key mapping for all operating systems:
KbName('UnifyKeyNames');

if ~exist('screenid', 'var')
    screenid = [];
end

if ~isempty(screenid)
    % 2nd argument exists: Numeric screenid or name of input calibfile?
    if isnumeric(screenid)
        % Screen id: Assign it.
        if length(screenid)~=1
            error('The "screenid" parameter must be a single integral screen number.');
        end
        
        if ~any(ismember(screenid, Screen('Screens')))
            error('Invalid "screenid" provided - no such screen.');
        end
        
        % Ok, seems to be a valid screenid. Null-Out the calibinfilename.
        calibinfilename = [];
        scal = [];
    else
        % Not a number, so it must be the filename of a calibration input
        % file:
        if ~ischar(screenid)
            error('You must provide the namestring of a valid calibration input file for parameter "calibinfilename"!');
        end
        
        % Namestring provided: Assign it, and Null-Out screenid.
        calibinfilename = screenid;
        screenid = [];
        
        % Check for such a file:
        if ~exist(calibinfilename, 'file')
            error(sprintf('No such calibration input file %s exists!', calibinfilename));
        end
        
        % Ok, valid file provided. Read it in:
        try
            load(calibinfilename);
        catch
            error(sprintf('Failed to load input calibration file %s for some reason.', calibinfilename));
        end
        
        fprintf('\nLoaded initial calibration data from file %s.\n', calibinfilename);
        
        if ~exist('warptype', 'var')
            fprintf('This file is not yet in PTB-3 format. It will be upgraded during save operation.\n');
        else
            if ~strcmp(warptype, 'BVLDisplayList')
                error('The calibration file is not a BVL calibration file, but some other incompatible format! Check your filenames - Aborted.');
            end
        end
        
        if ~exist('scal', 'var')
            error('The calibration file is not a BVL calibration file, but some other incompatible format! Check your filenames - Aborted.');
        end
        
        % Ok, warptype may or may not be initialized, but that doesn't matter.
        % 'scal' should contain initial calibration data...
        
        % Assign screenid from scal:
        screenid = scal.screenNumber;
        
        if isfield(scal, 'stereoMode')
            % Assign stereo mode from scal:
            stereomode = scal.stereoMode;
        end
    end
else
    % No 'screenid' provided: No calibration input file and no screenid:
    calibinfilename = [];
    scal = [];
    screenid = [];
end

% Target screen selection, if not yet defined by cmd-line input or calib
% input file:
if isempty(screenid)
    % Determine the arrangement of monitors, if more than one monitor is
    % used:
    scrns = Screen('Screens');

    % Decide what monitor to calibrate. This won't be asked if there is only
    % one monitor
    if length(scrns)>1,
        screenid = str2double(input(['Screen number of display to calibrate? [' num2str(scrns) '] '],'s'));
        if ~ismember(screenid, scrns)
            error('No such screen - Provide a valid id please! Aborted.');
        end
    else
        screenid = scrns;
    end    
end

if ~exist('stereomode', 'var')
    stereomode = [];
end

if isempty(stereomode)
    stereomode = 0;
end

% At this point, screenid contains the final screenid for the screen to
% calibrate.
% Assign it to scal struct. This will create 'scal' if it doesn't exist
% yet, or override its screenid in some cases:
scal.screenNumber = screenid;

% Ditto for stereoMode:
scal.stereoMode = stereomode;

% Fetch screen resolution: This is just used to build default output file
% name if none provided:
[w, h] = Screen('WindowSize', screenid);

% Check if name for calibration result file is provided:
if ~exist('caliboutfilename', 'var')
    caliboutfilename = [];
end

if isempty(caliboutfilename)
    % Nope: Assign default name - Store in dedicated subfolder of users PTB
    % config dir, with a well defined name that also encodes the screenid
    % and resolution for which to calibrate:
    caliboutfilename = [ PsychtoolboxConfigDir('GeometryCalibration') 'BVLCalibdata' sprintf('_%i_%i_%i', screenid, w, h) '.mat'];
    fprintf('\nNo name for calibration file provided. Using default name and location...\n');
end

% Print name of calibfile and check for existence of file:
fprintf('Name of calibration result file: %s\n\n', caliboutfilename);
if exist(caliboutfilename, 'file')
    answer = input('This file already exists. Overwrite it [y/n]? ','s');
    if ~strcmp(lower(answer), 'y') %#ok<STCI>
        fprintf('\n\nCalibration aborted. Please choose a different name for calibration result file.\n\n');
        return;
    end
end

% No need for synctests here...
oldsynclevel = Screen('Preference', 'SkipSyncTests', 2);
% Reduce PTB output to errors only:
oldverbosity = Screen('Preference', 'Verbosity', 1);

% Open the onscreen window for calibration on screen screenid:
[scal.windowPtr, screenrect] = Screen('OpenWindow', scal.screenNumber, 0, [], [], [], scal.stereoMode);

if isfield(scal, 'rect')
    % rect defined. Check for match with current screens rect:
    if sum(screenrect - scal.rect)~=0,
        Screen('CloseAll');
        error('The resolution of the monitor does not match the resolution stored in spatial calibration input file! Aborted.')
    end
else
    % rect not yet defined. Assign it:
    scal.rect = screenrect;
end

% Select a cross hair cursor:
ShowCursor('CrossHair');

% 'scal' has at least fields screenNumber, windowPtr and screen rect...

success = 0;

if ~exist('referenceImage', 'var')
    referenceImage = [];
end

if ~isempty(referenceImage)
    try
        refImg = imread(referenceImage);
    catch
        Screen('CloseAll');
        error(sprintf('Reference image file %s failed to load - No such file or permission problems?! Aborted.', referenceImage));
    end
    
    % Build reference texture and assign it to scal:
    scal.refTex = Screen('MakeTexture', scal.windowPtr, refImg);
end

% Define type of mapping for this calibration method:
% This is used in the CreateDisplayWarp() routine when parsing the
% calibration file to detect the type of undistortion method to use, ie.
% how to interpret the data in the calibration file to setup the
% calibration.
warptype = 'BVLDisplayList'; %#ok<NASGU>

% All following calibration code is try-catch protected...
try
    % Input file with starting point for calibration provided?
    if isempty(calibinfilename)
        % No. Create default calibration grid from scratch:
        
        % These are the initial settings for horizontal and vertical number
        % of evenly spaced calibration dots as used for the Banks Labs
        % haploscope:
        scal.NxdotsG = 37;
        scal.NydotsG = 27;

        if ~exist('xnum', 'var')
            xnum = [];
        end
        
        if ~exist('ynum', 'var')
            ynum = [];
        end
        
        % Any override values provided?
        if ~isempty(xnum)
            scal.NxdotsG = max(xnum, 1);
        end
        
        if ~isempty(ynum)
            scal.NydotsG = max(ynum, 1);
        end

        % Enforce odd number of points in new calibration interface
        if ~mod(scal.NxdotsG, 2)
            scal.NxdotsG = scal.NxdotsG + 1;
        end
        if ~mod(scal.NydotsG, 2)
            scal.NydotsG = scal.NydotsG + 1;
        end

        % Create default grid from scratch with full grid of calibration
        % dots:
        [scal] = createnewcalibrationgrid(scal);

        % Interactive selection of the dots that will be fitted. The dots x and y coordinates
        % will be stored in scal.FITDOTLIST
        [scal] = bvlSelectFitPts(scal);
        
        % Save intermediate calibration variables to file 'caliboutfilename'. This
        % method should work on both, Matlab 6.x, 7.x, ... and GNU/Octave - create
        % files that are readable by all runtime environments:
        save(caliboutfilename, 'warptype', 'scal', '-mat', '-V6');
        
        % MK: Superseded by direct call to bvlSelectFitPts above: [scal] = selectcalibrationpoints(scal);
    else
        % Yes. 'scal' should contain a good starting point. Take 'scal' as
        % starting point for interactive (re-)calibration procedure:

        % MK: Nothing to do yet as codepath shared with if-then branch...
    end

    % At this point, the initial 'scal' calibration grid is ready for
    % actual interactive calibration...
    
    % 'scal' contains initial calib params and grid setup. Perform interactive
    % GUI based calibration, adjustment and fitting:    
    scal.isDONTSTOP = 0;

    % Try to disable keystroke input to Matlab, but just continue if it fails:
    try
        ListenChar(2);
    catch
    end

    scal = manual_calibrate(scal);

    % Try to reenable keystroke input to Matlab, but just continue if it fails:
    try
        ListenChar(0);
    catch
    end

    % MK: Superseded by call to manual_calibrate above: scal = fineadjust(scal);

    % Save intermediate calibration variables to file 'caliboutfilename'. This
    % method should work on both, Matlab 6.x, 7.x, ... and GNU/Octave - create
    % files that are readable by all runtime environments:
    save(caliboutfilename, 'warptype', 'scal', '-mat', '-V6');

    % Interactive calibration finished, 'scal' is ready for numeric
    % polynomial fit to create the final useable calibration:
    
    % Perform fit operation:
    [scal] = fit(scal);
    
    % Successfully finished!
    success = 1;
catch
    % Close screen, restore everything...
    
    % Select standard arrow cursor:
    ShowCursor('Arrow');

    sca;
    Screen('Preference', 'SkipSyncTests', oldsynclevel);
    Screen('Preference', 'Verbosity', oldverbosity);
    
    % Rethrow the error for debugging...
    psychrethrow(psychlasterror);

    % Try to reenable keystroke input to Matlab, but just continue if it fails:
    try
        ListenChar(0);
    catch
    end
    
    % We try to continue, so the intermediate state gets saved to file as a
    % "restore point".
end

% Select standard arrow cursor:
ShowCursor('Arrow');

% Close display:
Screen('CloseAll');

% Restore sync test settings:
Screen('Preference', 'SkipSyncTests', oldsynclevel);
Screen('Preference', 'Verbosity', oldverbosity);

% 'scal' contains the final results of calibration. Write it out to
% calibfile for later use by the runtime routines:

% Save all relevant calibration variables to file 'caliboutfilename'. This
% method should work on both, Matlab 6.x, 7.x, ... and GNU/Octave - create
% files that are readable by all runtime environments:
save(caliboutfilename, 'warptype', 'scal', '-mat', '-V6');

if success
    fprintf('Calibration finished :-)\n\n');
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
    fprintf('The ''GeometryCorrection'' call has three additional optional parameters:\n');
    fprintf('The first one is a ''debug'' flag: Set it to a non-zero value for diagnostic output at runtime.\n');
    fprintf('The 2nd and 3rd one are override values for the xLoomSize and yLoomSize of the "warp mesh".\n');
    fprintf('Higher values mean a finer calibration but more computation time. Defaults are 73 and 53.\n');
    fprintf('Example:\n');
    fprintf('PsychImaging(''AddTask'', ''LeftView'', ''GeometryCorrection'', ''%s'', 0, 150, 100);\n', caliboutfilename);
    fprintf('would choose no debug output (=0), but a xLoomSize of 150 instead of 73 and a yLoomSize of 100\n');
    fprintf('instead of the default 53.\n');
    fprintf('For most purpose, the default settings will be fine though.\n\n');
    
else
    fprintf('CALIBRATION ABORTED DUE TO ERROR!\n\n');
    fprintf('The current state of calibration has been saved to %s.\n', caliboutfilename);
    fprintf('After fixing the issue, you can restart, using that file as "calibinfilename",\n');
    fprintf('so you don''t lose all of your work.\n\n');
end

% Done.
return;
