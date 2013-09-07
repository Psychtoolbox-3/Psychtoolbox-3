function DisplayUndistortionBezier(caliboutfilename, xnum, ynum, subdivision, imagename, screenid, stereomode, winrect, calibinfilename, refimagename)
% DisplayUndistortionBezier([caliboutfilename] [, xnum=2][, ynum=2][, subdivision=100][, imagename=default][, screenid=max][, stereomode=0][, winrect=[]][, calibinfilename][, refimagename])
%
% Geometric display calibration procedure for geometric undistortion of
% distorted displays. Needs graphics hardware with support for PTB imaging
% pipeline.
%
% IMPORTANT: While this routine is easy to use, it is also limited! The
% Banks Vision Lab (University of California, Berkeley) has contributed a
% much more powerful and flexible calibration procedure.
%
% Read "help DisplayUndistortionBVL" for help and info on that one.
%
%
% Many display devices, e.g., video beamers and most CRT displays, cause
% some amount of spatial distortion to your visual stimuli during display.
% Psychtoolbox can "undistort" your visual stimuli for you: At stimulus
% onset time, PTB applies a geometric warping transformation to your
% stimulus which is meant to counteract or cancel out the geometric
% distortion caused by your display device. If both PTB's warp transform
% and the distortion transform of the display match, your stimulus will
% show up undistorted on the display device.
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
% interactive procedure: It applies some warp transformation, shows the
% result on your display, and you can change and tweak the transformation
% until it fits your needs, ie. it nicely undistorts your display. Then the
% corresponding calibration file is saved for later use with that display.
%
% As the name suggests, DisplayUndistortionBezier supports a continous
% mapping (x', y') = Beziersurface(x, y) from input pixel locations (x,y)
% in your stimulus image to output locations (x', y') on your display. This
% mapping is defined by a Bezier surface (see Chapter 12 "Evaluators and
% NURBS" or any other textbook about NURBS for a description of their
% mathematical properties).
%
% The shape of the Bezier surface is defined by the location of a couple of
% control points, which you can move around during the calibration
% procedure to modify the shape of the surface. You can select the number
% of control points to use: More control points (=degrees of freedom) allow
% for more flexibility and finer control, but make the calibration procedure
% more tedious for you. Their number doesn't affect computation time when
% your experiment script is running.
%
% How to use:
% 
% 1. Start the script, providing all parameters that you don't want to have
% at default settings (all parameters have defaults):
%
% 'caliboutfilename' Name of the file to which calibration results should
% be stored. Defaults to 'BezierCalibdata.mat' in your current working directory.
%
% 'xnum' and 'ynum' Number of horizontal and vertical control points to
% use. Higher numbers mean finer and more flexible control, but also a more
% tedious calibration for you. Their number doesn't affect runtime
% behaviour of your stimulus script. For simple trapezoid correction or
% translation/rotation, the minimum allowed number of control points xnum=2
% and ynum=2 is sufficient - this is the default. For most other purposes
% a 3 by 3 grid of control points xnum=3 and ynum=3 may suffice. However,
% you are not limited by any upper bound...
%
% 'subdivision' Number of vertical and horizontal subdivisions of the
% bezier surface - the grid resolution: Higher numbers mean higher accuracy
% but also higher computational overhead in your script. However, recent
% graphics hardware shouldn't have much problems handling reasonably sized
% meshes. Defaults to a 100 by 100 grid.
%
% 'imagename' Name of the image file for the test image to be rendered as
% an alternative to the mesh grid. We default to our standard 'konintjes'
% image.
%
% 'screenid' screen id of the target display for calibration.
% max(Screen('Screens')) by default.
%
% 'stereomode' Stereomode for which calibration should be applied: Defaults
% to 0 == Mono mode. [6 1] would mean: "Use stereomode 6 (Anaglyph stereo)
% and the right-eye view (1)".
%
% 'winrect' Size of the calibration window: Defaults to full-screen.
%
% 'calibinfilename' Defaults to none. If provided, results of a previous
% calibration are loaded from file 'calibinfilename' instead of starting
% from scratch. Useful for incremental calibration.
%
% 2. After startup, the script will display a grid onscreen, which
% represents the displayed area after calibration. Your job is to tweak,
% shift and bend that grid so that it looks as flat and rectilinear as
% possible on your display from the viewpoint of your subject. The grid has
% green control points placed at regular intervals. These are the tweakable
% points that you can move: Moving these control points will bend the
% calibration grid in a smooth fashion, as if the grid would be attached to
% the points with some springs.
%
% Mouse operation:
%
% To select a control point, just move the mouse pointer close to it, then
% press a mouse key. The selected control point will change color to yellow
% and a yellow line will connect its current position to its original
% position in the uncalibrated display.
%
% You can move the point by moving the mouse pointer while keeping the
% mouse key pressed.
%
% 
% Keys and their meaning:
%
% You can also move the currently selected control point via the Cursor
% keys, at slow speed, or at a faster speed when holding down the shift
% key.
%
% Press the 'space' key to toggle between grid display and display of the
% test image 'imagename'. A good way to test the calibration would be to
% load a screenshot of one of your typical stimuli as image file
% 'imagename'.
%
% You finish the calibration and write it into a calibration file by
% pressing the ESCape key. This will end the calibration script.

% History:
% 07/16/07 Initial (incomplete) version. (MK)

% Running on PTB-3? Abort otherwise:
PsychDefaultSetup(1);

% No need for synctests here...
oldsynclevel = Screen('Preference', 'SkipSyncTests', 1);

KbName('UnifyKeyNames');
UpArrow = KbName('UpArrow');
DownArrow = KbName('DownArrow');
LeftArrow = KbName('LeftArrow');
RightArrow = KbName('RightArrow');
esc = KbName('ESCAPE');
shift = KbName('RightShift');
space = KbName('space');

if nargin < 1 || isempty(caliboutfilename)
    caliboutfilename = 'BezierCalibdata.mat';
    fprintf('Warning: No name for calibration file provided. Using default name...\n');
end
fprintf('Name of calibration result file: %s\n', caliboutfilename);
if exist(caliboutfilename, 'file')
    answer = input('This file already exists. Overwrite it [y/n]?','s');
    if ~strcmp(lower(answer), 'y') %#ok<STCI>
        error('Calibration aborted. Please choose a different name for calibration result file.');
    end
end

if nargin < 2 || isempty(xnum)
    xnum = 2;
end

if xnum < 2
    fprintf('xnum set to an invalid value (smaller than 2): Reset to 2...\n');
    xnum = 2;
end

if nargin < 3 || isempty(ynum)
    ynum = 2;
end

if ynum < 2
    fprintf('ynum set to an invalid value (smaller than 2): Reset to 2...\n');
    ynum = 2;
end

if nargin < 4 || isempty(subdivision)
    subdivision = 100;
end

fprintf('Will use %i subdivisions of calibration warp mesh for undistortion.\n', subdivision);


if nargin < 5 || isempty(imagename)
    imagename = [PsychtoolboxRoot 'PsychDemos/konijntjes1024x768.jpg'];
end

if nargin < 6 || isempty(screenid)
    screenid = max(Screen('Screens'));
end

if nargin < 7 || isempty(stereomode)
    viewid = 0;
    stereomode = 0;
else
    viewid = stereomode(2);
    stereomode = stereomode(1);
end

if nargin < 8 || isempty(winrect)
    winrect = [];
end

if nargin < 9
    calibinfilename = [];
end

if nargin < 10
    refimagename = [];
end

if ~isempty(refimagename)
    refimg = imread(refimagename);    
end

if nargin >= 8 && ~isempty(winrect) && strcmpi(winrect, 'REFIMAGE')
    if ~isempty(refimagename)
        winrect = [0, 0, size(refimg, 2), size(refimg, 1)];
    else
        winrect = [];
    end
end

InitializeMatlabOpenGL([], [], 1);
win = Screen('OpenWindow', screenid, 0, winrect, [], [], stereomode, [], mor(kPsychNeedFastBackingStore, kPsychNeedOutputConversion));

if ~isempty(refimagename)
    reftex = Screen('MakeTexture', win, refimg);
else
    reftex = [];
end

% Allocate display list handle and build initial warpstruct:
gld = glGenLists(1);
warpstruct.glsl = [];
warpstruct.gld = gld;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'FinalFormatting', 'GeometryCorrection', warpstruct);
PsychImaging('FinalizeConfiguration');
PsychImaging('PostConfiguration', win);

img = imread(imagename);
tex = Screen('MakeTexture', win, img);

Screen('Flip', win);
Screen('SelectStereoDrawBuffer', win, viewid);

% Determine width and height of window:
[w, h] = Screen('WindowSize', win);

% Read parameters from previous calibration file?
if ~isempty(calibinfilename)
    % Old file provided: Setup from that file:
    calib = load(calibinfilename);
    if ~strcmp(calib.warptype, 'BezierDisplayList')
        Screen('CloseAll');
        error('Provided input calibration file %s does not describe a calibration created with this routine!', calibinfilename);
    end
    
    subdivision = calib.subdivision;
    frompts = calib.frompts;
    topts = calib.topts;
else
    % No old file: Initialize based on call arguments:

    % Setup initial mapping table for texture coordinates (source image control
    % points):
    frompts = zeros(2, 2, 2);
    frompts(:, 1, 1) = [0 h];
    frompts(:, 1, 2) = [w h];
    frompts(:, 2, 2) = [w 0];
    frompts(:, 2, 1) = [0 0];

    % Setup initial mapping table for vertex coordinates (target image control
    % points). We start with a uniform rectilinear spacing of points:
    topts = zeros(3, ynum, xnum);
    dx = w / (xnum - 1);
    dy = h / (ynum - 1);

    for y=1:ynum
        for x=1:xnum
            topts(1, y, x) = ((x-1) * dx);
            topts(2, y, x) = ((y-1) * dy);
        end
    end
end

% Setup a 2D parametric grid with 'subdivision' subdivisions:
glMapGrid2d(subdivision, 0, 1, subdivision, 0, 1);

% Enable Bezier evaluators:
glEnable(GL.MAP2_VERTEX_3);
glEnable(GL.MAP2_TEXTURE_COORD_2);

% Establish mapping for texture coordinates:
glMap2d(GL.MAP2_TEXTURE_COORD_2, 0, 1, 2, size(frompts,2), 0, 1, 2*size(frompts,2), size(frompts,3), frompts);

% Make sure all mouse buttons are released:
oldbuttons = 1;
while oldbuttons
    [xm ym oldbuttons] = GetMouse(win);
    oldbuttons = any(oldbuttons);
end

maxx = -1;
maxy = -1;
oldtopts = topts;
applycalib = 0;

% Prevent keypresses from spilling into Matlab window:
ListenChar(2);

% Calibration loop: Runs until keypress:
while 1
    % Setup mapping based on current control point matrix for destination
    % points:
    glMap2d(GL.MAP2_VERTEX_3, 0, 1, 3, size(topts,2), 0, 1, 3*size(topts,2), size(topts,3), topts);

    % Drawing color for calibration grid is red:
    glColor3f(1, 0, 0);

    if applycalib
        % Build a display list that corresponds to the current calibration:
        glNewList(gld, GL.COMPILE);
        
        % Compute the mesh based on current mappings:
        glEvalMesh2(GL.FILL, 0, subdivision, 0, subdivision);
        
        % List ready - and already updated in the imaging pipeline:
        glEndList;
        
        % Draw test image to backbuffer:
        Screen('FrameRect', win, [255 0 0], [1 1 w-1 h-1]);
        Screen('DrawTexture', win, tex);
        
        % Enable pipeline: It will use the computed calibration:
        Screen('HookFunction', win, 'Enable', 'FinalOutputFormattingBlit');
        
        % The next Screen('Flip') will show the testimage with calibration
        % applied...
    else
        % Visualization of the calibration grid requested.

        % Draw backdrop image, if there is one:
        if ~isempty(reftex)
            Screen('DrawTexture', win, reftex, [], Screen('Rect', reftex));

            % Reset drawing color for calibration grid to red:
            glColor3f(1, 0, 0);
        end

        % Draw calibration grid:
        glEvalMesh2(GL.LINE, 0, subdivision, 0, subdivision);
        
        % Disable pipeline so the mesh gets shown one-to-one onscreen:
        Screen('HookFunction', win, 'Disable', 'FinalOutputFormattingBlit');
    end
    
    % Draw control points in green:
    for y=1:size(topts, 2)
        for x=1:size(topts, 3)
            Screen('DrawDots', win, topts(1:2, y, x), 10, [0 255 0]);
        end
    end
    
    % Draw currently selected control point and its connecting line to
    % initial position in yellow:
    if maxx~=-1
        Screen('DrawDots', win, topts(1:2, maxy, maxx), 10, [255 255 0]);
        Screen('DrawLine', win, [255 255 0], topts(1, maxy, maxx), topts(2, maxy, maxx), oldtopts(1, maxy, maxx), oldtopts(2, maxy, maxx));
    end
    
    % Show it:
    Screen('Flip', win);
    
    % Query mouse:
    [xm ym buttons] = GetMouse(win);
    buttons = any(buttons);
    
    if buttons
        % Mouse button pressed. First press?
        if oldbuttons
            % Nope. This is an update of the control point location:
            topts(1, maxy, maxx) = xm;
            topts(2, maxy, maxx) = ym;
        else
            % Yes. Find closest controlpoint and select it as active one:
            maxdist = inf;
            maxx = -1;
            maxy = -1;
            for y=1:size(topts, 2)
                for x=1:size(topts, 3)
                    dx = xm - topts(1, y, x);
                    dy = ym - topts(2, y, x);
                    d = sqrt(dx*dx + dy*dy);
                    if d < maxdist
                        maxdist = d;
                        maxx = x;
                        maxy = y;
                    end
                end
            end            
        end
    end

    % Store old buttons state:
    oldbuttons = buttons;

    % Check keyboard:
    [isdown secs keycode] = KbCheck;
    if isdown
        if maxx~=-1
            % Holding down shift increases the stepwidth of keyboard based
            % control point movements:
            if keycode(shift)
                delta = 1;
            else
                delta = 0.1;
            end
            
            % Movement of current selected controlpoint via cursor keys:
            if keycode(UpArrow)
                topts(2, maxy, maxx) = topts(2, maxy, maxx) - delta;
            end
            if keycode(DownArrow)
                topts(2, maxy, maxx) = topts(2, maxy, maxx) + delta;
            end
            if keycode(LeftArrow)
                topts(1, maxy, maxx) = topts(1, maxy, maxx) - delta;
            end
            if keycode(RightArrow)
                topts(1, maxy, maxx) = topts(1, maxy, maxx) + delta;
            end
        end
        
        % Space toggles between grid visualization and visualization of
        % applied calibration:
        if keycode(space)
            applycalib = 1 - applycalib;
            while KbCheck; end;
        end
        
        % ESCape key marks end of calibration procedure:
        if keycode(esc)
            break;
        end
    end
end

% End of calibration.

% Reenable Matlabs keyboard handling:
ListenChar(0);

% Reenable sync tests:
Screen('Preference', 'SkipSyncTests', oldsynclevel);

% Disable mesh evaluators:
glDisable(GL.MAP2_VERTEX_3);
glDisable(GL.MAP2_TEXTURE_COORD_2);

% Writeout of calibration/undistortion matrices:

% Define type of mapping for this calibration method:
warptype = 'BezierDisplayList'; %#ok<NASGU>

% Save relevant calibration variables to file 'caliboutfilename':
save(caliboutfilename, 'warptype', 'subdivision', 'frompts', 'topts', '-mat', '-V6');

% Close Display:
Screen('CloseAll');

% Done.
return;
