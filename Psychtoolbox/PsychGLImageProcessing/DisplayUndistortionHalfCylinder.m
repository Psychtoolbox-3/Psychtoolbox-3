function DisplayUndistortionHalfCylinder(calibfilename, screenid)
% Create geometric display calibration file for projection onto a cylinder or sphere.
%
% CAUTION: This function is not fully implemented and tested! It may or
% may not work for your purpose and may give totally wrong results.
%
% Usage:
%
% DisplayUndistortionHalfCylinder([calibfilename][, screenid]);
%
% Prepare display calibration file for projection of a flat screen image to
% a half-cylindrical projection surface.
%
% The script displays a calibration grid on screen 'screenid' (or the
% auto-selected screen if this parameter is omitted). Use the keys
% mentioned below to modify the grid until it is optimally aligned with
% your cylindrical projection surface. Save the final calibration settings
% to file 'calibfilename' (or an auto-selected default filename, if
% parameter is omitted) by pressing the 's' key. Exit the procedure via
% pressing the ESCape key.
%
% The saved calibration file can be later loaded and applied in a script
% via the command sequence that this script will diplay when it is
% finished.
%
%
% Keys and their meaning:
% ESCape = Quit script. Unsaved changes will be lost!
% s      = Save current settings in calibration file.
% r      = Change rotation angle.
% t      = Change translation of output image.
% o      = Change size of output image.
% i      = Change size of input image.
% w      = Width of the projection area at a radius distance R in cm for
%          sphere curvature correction.
% b      = Radius of the spherical projection screen for sphere curvature correction.
%
% -> Type "edit SphereProjectionShader.frag.txt" and look at the warp
% shader code for a full understanding of the roff and rpow coefficients
% for spherical remapping. Also see Psychtoolbox forum message 11660.
%
%
% Arrow keys change the parameter that is selected via 'r', 't', 'o' or 'i'
% key.
%

% History:
% 13.4.2009  mk Written.
% 19.12.2010 mk Updated with roff, rpow parameter setting for spherical
%               projections as suggested by Ingmar Schneider.
% 25.8.2011  mk Adapt code for sphere projection undistortion to new
%               convention of Ingmar Schneider's shader code.

global GL;

PsychDefaultSetup(1);

correctiontype = input('Cylinder projection (c) or sphere projection (s)? ', 's');
if ~ismember(correctiontype, {'c','s'})
    error('Invalid correction type provided! Must be ''s'' or ''c''.'); 
end

KbName('UnifyKeyNames');
UpArrow = KbName('UpArrow');
DownArrow = KbName('DownArrow');
LeftArrow = KbName('LeftArrow');
RightArrow = KbName('RightArrow');
esckey = KbName('ESCAPE');
rkey = KbName('r');
tkey = KbName('t');
okey = KbName('o');
ikey = KbName('i');
skey = KbName('s');
Wflatkey = KbName('w');
Rkey = KbName('b');

if ~exist('screenid', 'var') || isempty(screenid)
    screenid=max(Screen('Screens'));
end

if ~exist('calibfilename', 'var')
    calibfilename = [];
end

if isempty(calibfilename)
    % Nope: Assign default name - Store in dedicated subfolder of users PTB
    % config dir, with a well defined name that also encodes the screenid
    % and resolution for which to calibrate:
    [w, h] = Screen('WindowSize', screenid);
    
    if correctiontype == 'c'
        defname = 'HalfCylinderCalibdata';
    else
        defname = 'SphereCalibdata';
    end
    
    calibfilename = [ PsychtoolboxConfigDir('GeometryCalibration') defname sprintf('_%i_%i_%i', screenid, w, h) '.mat'];
    fprintf('\nNo name for calibration file provided. Using default name and location...\n');
end

% Print name of calibfile and check for existence of file:
fprintf('Name of calibration result file: %s\n\n', calibfilename);
if exist(calibfilename, 'file')
    answer = input('This file already exists. Overwrite it [y/n]? ','s');
    if ~strcmp(lower(answer), 'y') %#ok<STCI>
        fprintf('\n\nCalibration aborted. Please choose a different name for calibration result file.\n\n');
        return;
    end
end

% Create initial default cylinder undistortion calibration struct:
if correctiontype == 'c'
    calib.warptype = 'HalfCylinderProjection';
else
    calib.warptype = 'SphereProjection';
end

calib.rotationAngle = [];
calib.inSize = [];
calib.inOffset = [];
calib.outOffset = [];
calib.outSize = [];
calib.Wflat = [];
calib.R = [];

firsttime = 1;

try
    ListenChar(2);

    % Open onscreen window, with this calibration applied:

    % Prepare pipeline for configuration. This marks the start of a list of
    % requirements/tasks to be met/executed in the pipeline:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'AllViews', 'GeometryCorrection', calib);
    win =PsychImaging('OpenWindow', screenid, 0);
    Screen('TextSize', win, 16);

    % Alpha blending for nice anti-aliased calibration dot grid:
    Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    % Assign reasonable defaults for all calibration parameters:
    [winWidth, winHeight] = Screen('WindowSize', win);
    calib.rotationAngle = 0;
    calib.inSize = [winWidth, winHeight];
    calib.inOffset = [0, 0];
    calib.outOffset = [0, 0];
    calib.outSize = [winWidth, winHeight];
    
    % These screen width, and sphere radius settings are reasonable
    % defaults, used by Ingmar Schneider for a tested spherical display
    % setup:
    calib.Wflat = 44.5;
    calib.R = 32.0392;
    
    % Query generated warpstruct:
    warpstruct = CreateDisplayWarp('Query');

    % Build dotgrid for visualization of current structure:
    [y,x] = meshgrid(0:10:winHeight, 0:10:winWidth);
    xy = [x(:)' ; y(:)'];

    xylines = zeros(2,0);
    li = 0;
    for i = 1:length(y)
        li = li + 10;
        xylines = [ xylines, [1, li]', [winWidth-1, li]' ]; %#ok<AGROW>
    end

    for i = 1:length(x)
        li = x(i);
        xylines = [ xylines, [li, 0]', [li, winHeight]' ]; %#ok<AGROW>
    end

    mode = 1;

    % Visualization and recalibration loop:
    while 1
        % Draw calibration grid:
        Screen('DrawLines', win, xylines, 1, [0 255 0], [],1);
        Screen('DrawDots', win, xy, 4, [255 255 0], [], 1);
        % Show updated calibration:
        Screen('Flip',win);

        % Keypress?
        [down, secs, keyCode] = KbCheck;
        if down || firsttime
            if keyCode(skey)
                % Save current calibration:
                warptype = calib.warptype; %#ok<NASGU>
                rotationAngle = calib.rotationAngle; %#ok<NASGU>
                inSize = calib.inSize; %#ok<NASGU>
                inOffset = calib.inOffset; %#ok<NASGU>
                outOffset = calib.outOffset; %#ok<NASGU>
                outSize = calib.outSize; %#ok<NASGU>
                Wflat = calib.Wflat; %#ok<NASGU>
                R = calib.R; %#ok<NASGU>
                
                % Save all relevant calibration variables to file 'calibfilename'. This
                % method should work on both, Matlab 6.x, 7.x, ... and GNU/Octave - create
                % files that are readable by all runtime environments:
                save(calibfilename, 'warptype', 'rotationAngle', 'inSize', 'inOffset', 'outOffset', 'outSize', 'Wflat', 'R', '-mat', '-V6');

                Beeper;
                continue;
            end

            if keyCode(esckey)
                break;
            end

            if keyCode(rkey)
                mode = 1;
            end
            if keyCode(tkey)
                mode = 2;
            end
            if keyCode(okey)
                mode = 3;
            end
            if keyCode(ikey)
                mode = 4;
            end
            if keyCode(Wflatkey)
                mode = 5;
            end
            if keyCode(Rkey)
                mode = 6;
            end

            dx = 0;
            if keyCode(LeftArrow)
                dx = -1;
            end
            if keyCode(RightArrow)
                dx = +1;
            end

            dy = 0;
            if keyCode(UpArrow)
                dy = -1;
            end
            if keyCode(DownArrow)
                dy = +1;
            end

            switch(mode)
                case 1,
                    calib.rotationAngle = calib.rotationAngle + dx / 10;
                case 2,
                    calib.outOffset = calib.outOffset + [dx, dy];
                case 3,
                    calib.outSize = calib.outSize  + [dx, dy];
                case 4,
                    calib.inSize = calib.inSize  + [dx, dy];
                case 5,
                    old = calib.Wflat;
                    Screen('HookFunction', win, 'Disable', 'StereoLeftCompositingBlit');
                    calib.Wflat = GetEchoNumber(win, sprintf('Enter width of the projection area in cm at the radius distance R in cm for sphere proj. [%f]: ', calib.Wflat), 30, 30, 255, []);
                    Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                    if isempty(calib.Wflat)
                        calib.Wflat = old;
                    end
                    mode = -1;
                case 6,
                    old = calib.R;
                    Screen('HookFunction', win, 'Disable', 'StereoLeftCompositingBlit');
                    calib.R = GetEchoNumber(win, sprintf('Enter radius of the used spherical screen in cm for sphere proj. [%f]: ', calib.R), 30, 30, 255, []);
                    Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');
                    if isempty(calib.R)
                        calib.R = old;
                    end
                    mode = -1;
                otherwise
                    continue
            end

            % Build the unwarp mesh display list within the OpenGL context of
            % Psychtoolbox:
            Screen('BeginOpenGL', win, 1);

            % Build a display list that corresponds to the current calibration:
            glNewList(warpstruct.gld, GL.COMPILE);

            % "Draw" the warp-mesh once, so it gets recorded in the display
            % list:

            % No color gain correction:
            glColor4f(1,1,1,1);

            glTranslatef(calib.outOffset(1), calib.outOffset(2), 0);

            % Apply some rotation correction for misaligned displays:
            glTranslatef(calib.outSize(1)/2, calib.outSize(2)/2, 0);
            glRotatef(calib.rotationAngle, 0.0, 0.0, 1.0);
            glTranslatef(-calib.outSize(1)/2, -calib.outSize(2)/2, 0);

            % Draw a single default quad:
            glBegin(GL.QUADS)
            glTexCoord2f(0,calib.outSize(2));
            glVertex2f(0,0);

            glTexCoord2f(calib.outSize(1),calib.outSize(2));
            glVertex2f(calib.outSize(1),0);

            glTexCoord2f(calib.outSize(1),0);
            glVertex2f(calib.outSize(1),calib.outSize(2));

            glTexCoord2f(0,0);
            glVertex2f(0,calib.outSize(2));
            glEnd;

            % List ready - and already updated in the imaging pipeline:
            glEndList;

            Screen('EndOpenGL', win);

            % Update shader settings:
            glUseProgram(warpstruct.glsl);
            glUniform2f(glGetUniformLocation(warpstruct.glsl, 'inSize'), calib.inSize(1), calib.inSize(2));
            glUniform2f(glGetUniformLocation(warpstruct.glsl, 'inOffset'), calib.inOffset(1), calib.inOffset(2));
            glUniform2f(glGetUniformLocation(warpstruct.glsl, 'outSize'), calib.outSize(1), calib.outSize(2));
            if strcmpi(calib.warptype, 'SphereProjection')
                % Additional parameters for sphere projection:
                glUniform1f(glGetUniformLocation(warpstruct.glsl, 'Wflat'), calib.Wflat);
                glUniform1f(glGetUniformLocation(warpstruct.glsl, 'R'), calib.R);
            end
            glUseProgram(0);

            % Play tricks to get rid of visual artifacts:
            Screen('HookFunction', win, 'Disable', 'StereoLeftCompositingBlit');
            Screen('FillRect', win, 0);
            Screen('Flip', win);
            Screen('Flip', win);
            Screen('HookFunction', win, 'Enable', 'StereoLeftCompositingBlit');

            firsttime = 0;
        else
            WaitSecs(0.01);
        end
        % Next iteration.
    end

    % Done.
    Screen('CloseAll');
    ListenChar(0);

    fprintf('Calibration finished :-)\n\n');
    fprintf('You can apply the calibration in your experiment script by replacing your \n')
    fprintf('win = Screen(''OpenWindow'', ...); command by the following sequence of \n');
    fprintf('commands:\n\n');
    fprintf('PsychImaging(''PrepareConfiguration'');\n');
    fprintf('PsychImaging(''AddTask'', ''LeftView'', ''GeometryCorrection'', ''%s'');\n', calibfilename);
    fprintf('win = PsychImaging(''OpenWindow'', ...);\n\n');
    fprintf('This would apply the calibration to the left-eye display of a stereo setup.\n');
    fprintf('Additional options would be ''RightView'' for the right-eye display of a stereo setup,\n');
    fprintf('or ''AllViews'' for both views of a stereo setup or the single display of a mono\n');
    fprintf('setup.\n\n');

    return;
catch
    Screen('CloseAll');
    ListenChar(0);
    psychrethrow(psychlasterror);
end
