function Kinect3DDemo(stereomode)
% Kinect3DDemo - Capture and display video and depths data from a Kinect box.
%
% Usage:
%
% Kinect3DDemo([stereomode=0])
%
% This connects to a Microsoft Kinect device on the USB bus, then captures
% and displays video and depths data delivered by the Kinect.
%
% This is an early prototype.
%
% Control keys and their meaning:
%
% 'a' == Zoom out by moving object away from viewer.
% 'z' == Zoom in by moving object close to viewer.
% Left- and Right cursor keys == Rotate object around vertical axis.
% Up- and Down cursor keys == Rotate object around horizontal axis.
% SPACE = Toggle live 3D feed.
% 'm' = Toggle mesh display vs. point-cloud display.
% ESCAPE == Quit demo.
%
% Options:
%
% stereomode = n. For n > 0 this activates stereoscopic rendering - The shape is
% rendered from two slightly different viewpoints and one of Psychtoolbox's
% built-in stereo display algorithms is used to present the 3D stimulus. This
% is very preliminary so it doesn't work that well yet.
%
% See help PsychKinect, help PsychKinectCore and help InstallKinect for
% further info.
%

% History:
% 5.12.2010  mk  Written.

global win;

% Is the script running in OpenGL Psychtoolbox?
AssertOpenGL;

%PsychDebugWindowConfiguration;
dst1 = [0, 0, 640, 480];
dst2 = [650, 0, 650+640, 480];

if nargin < 1 || isempty(stereomode)
    stereomode = 0;
end;

% Response keys: Mapping of keycodes to keynames.
KbName('UnifyKeyNames');
closer = KbName('a');
farther = KbName('z');
quitkey = KbName('ESCAPE');
rotateleft = KbName('LeftArrow');
rotateright = KbName('RightArrow');
rotateup = KbName('UpArrow');
rotatedown = KbName('DownArrow');
toggleCapture = KbName('space');
toggleMesh = KbName('m');

% Find the screen to use for display:
screenid=max(Screen('Screens'));

% Disable Synctests for this simple demo:
Screen('Preference','SkipSyncTests',1);

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab wrapper. We need to do this before the first call
% to any OpenGL function:
%InitializeMatlabOpenGL(0,1);
InitializeMatlabOpenGL(0,0);

% Open a double-buffered full-screen window: Everything is left at default
% settings, except stereomode:
rect = [];
%rect = [0 0 1300 500];

if stereomode > 0
    PsychImaging('PrepareConfiguration');
    [win , winRect] = PsychImaging('OpenWindow', screenid, 0, rect, [], [], stereomode);
else
    [win , winRect] = Screen('OpenWindow', screenid, 0, rect, [], [], stereomode);
end

if stereomode > 5 && stereomode < 10
    SetAnaglyphStereoParameters('OptimizedColorAnaglyphMode', win);
end

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Get the aspect ratio of the screen, we need to correct for non-square
% pixels if we want undistorted displays of 3D objects:
ar=winRect(4)/winRect(3);

% Enable proper occlusion handling via depth tests:
glEnable(GL.DEPTH_TEST);

% Set projection matrix: This defines a perspective projection,
% corresponding to the model of a pin-hole camera - which is a good
% approximation of the human eye and of standard real world cameras --
% well, the best aproximation one can do with 3 lines of code ;-)
glMatrixMode(GL.PROJECTION);
glLoadIdentity;

% Field of view is +/- 20 degrees from line of sight. Objects close than
% 0.1 distance units or farther away than 200 distance units get clipped
% away, aspect ratio is adapted to the monitors aspect ratio:
gluPerspective(20.0,1/ar,0.1,10000.0);

% Setup modelview matrix: This defines the position, orientation and
% looking direction of the virtual camera:
glMatrixMode(GL.MODELVIEW);
glLoadIdentity;

% Set size of points for drawing of reference dots
glPointSize(3.0);

% Setup initial z-distance of objects:
zz = 50.0;
ang = 0.0;      % Initial rotation angle
tilt = 0.0;

% Half eye separation in length units for quick and dirty stereoscopic
% rendering. Our way of stereo is not correct, but it makes for a
% nice demo. Figuring out proper values is not too difficult, but
% left as an exercise to the reader.
eye_halfdist=3;
if stereomode == 0
    eye_halfdist=0;
end

% Finish OpenGL setup and check for OpenGL errors:
Screen('EndOpenGL', win);

% Retrieve duration of a single monitor flip interval: Needed for smooth
% animation.
ifi = Screen('GetFlipInterval', win);

kinect = PsychKinect('Open');

%PsychKinect('ApplyCalibrationFile', kinect, '/home/kleinerm/Downloads/kinect_calibration.yml');

% Enable this to test kinects motor:
if 0
    PsychKinect('SetAngle', kinect, -30);
    WaitSecs('YieldSecs', 0.5);
    for angle = -30:30
        PsychKinect('SetAngle', kinect, angle);
        WaitSecs('YieldSecs', 0.50);
    end

    % Bring kinect into upright level position:
    PsychKinect('SetAngle', kinect, 0);
end

% Start kinects capture:
PsychKinect('Start', kinect);

% Initially sync us to the VBL:
vbl=Screen('Flip', win);

% Some stats...
tstart=vbl;
framecount = 0;
waitframes = 1;

doMesh = 0;
doCapture = 1;
oldKeyIsDown = KbCheck;
tex1 = [];
tex2 = [];
mymesh = [];

% Animation loop: Run until key press or one minute has elapsed...
t = GetSecs;
while ((GetSecs - t) < 600)
    if doCapture
        % Release old mesh, if any:
        if 1 && ~isempty(mymesh)
            PsychKinect('DeleteObject', win, mymesh);
            mymesh = [];
        end

        % Try to get hold of a new frame from the device, block if none available:
        [rc, cts] = PsychKinect('GrabFrame', kinect, 1);
        if rc > 0
            % Output video image from rgb camera if enabled:
            if 0
                [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 1, 1);
                if width > 0 && height > 0
                    tex2 = Screen('SetOpenGLTextureFromMemPointer', win, tex2, imbuff, width, height, channels, 1, GL.TEXTURE_RECTANGLE_EXT);
                end
            else
                if ~isempty(tex2)
                    Screen('Close', tex2);
                end
                tex2 = [];
            end

            % Create 3D mesh object with 3D encoded scene as seen by kinect:
            mymesh = PsychKinect('CreateObject', win, kinect);

            % Done fetching data for this kinect frame: Release it.
            PsychKinect('ReleaseFrame', kinect);

            % No valid mesh? Retry above.
            if isempty(mymesh)
                continue;
            end
        else
            % No valid data from kinect? Retry above.
            continue;
        end
    end

    % Draw kinect color image if enabled:
    if 0 && ~isempty(mymesh.tex) && (mymesh.tex > 0)
        Screen('DrawTexture', win, mymesh.tex, [], dst1);
    end

    % Draw kinect depths image if enabled:
    if ~isempty(tex2) && (tex2 > 0)
        Screen('DrawTexture', win, tex2, [], dst2);
    end

    % Switch to OpenGL rendering for drawing of next frame:
    Screen('BeginOpenGL', win);

    % Left-eye cam is located at 3D position (-eye_halfdist,0,zz), points upright (0,1,0) and fixates
    % at the origin (0,0,0) of the worlds coordinate system:
    glLoadIdentity;
    gluLookAt(-eye_halfdist, 0, zz, 0, 0, 0, 0, 1, 0);

    % Draw into image buffer for left eye:
    Screen('EndOpenGL', win);
    Screen('SelectStereoDrawBuffer', win, 0);
    Screen('BeginOpenGL', win);

    % Clear out the depth-buffer for proper occlusion handling:
    glClear(GL.DEPTH_BUFFER_BIT);

    % Call our subfunction that does the actual drawing of the shape (see below):
    drawShape(mymesh, doMesh, ang, tilt);

    % Stereo rendering requested?
    if (stereomode > 0)
        % Yes! We need to render the same object again, just with a different
        % camera position, this time for the right eye:

        % Right-eye cam is located at 3D position (+eye_halfdist,0,zz), points upright (0,1,0) and fixates
        % at the origin (0,0,0) of the worlds coordinate system:
        glLoadIdentity;
        gluLookAt(+eye_halfdist, 0, zz, 0, 0, 0, 0, 1, 0);
        % Draw into image buffer for right eye:
        Screen('EndOpenGL', win);
        Screen('SelectStereoDrawBuffer', win, 1);
        Screen('BeginOpenGL', win);

        % Clear out the depth-buffer for proper occlusion handling:
        glClear(GL.DEPTH_BUFFER_BIT);

        % Call subfunction that does the actual drawing of the shape (see below):
        drawShape(mymesh, doMesh, ang, tilt);
    end;

    % Finish OpenGL rendering into Psychtoolbox - window and check for OpenGL errors.
    Screen('EndOpenGL', win);

    % Tell Psychtoolbox that drawing of this stim is finished, so it can optimize
    % drawing:
    Screen('DrawingFinished', win);

    % Now that all drawing commands are submitted, we can do the other stuff before
    % the Flip:

    % Check for keyboard press:
    [KeyIsDown, endrt, KeyCode] = KbCheck(-1);
    if KeyIsDown
        if ( KeyIsDown==1 && KeyCode(closer)==1 )
            zz=zz-0.25;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(farther)==1 )
            zz=zz+0.25;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotateright)==1 )
            ang=ang+1.0;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotateleft)==1 )
            ang=ang-1.0;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotatedown)==1 )
            tilt=tilt+1.0;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(rotateup)==1 )
            tilt=tilt-1.0;
            KeyIsDown=0;
        end
    end

    if KeyIsDown && ~oldKeyIsDown
        oldKeyIsDown = KeyIsDown;

        if ( KeyIsDown==1 && KeyCode(toggleCapture)==1 )
            doCapture = ~doCapture;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(toggleMesh)==1 )
            doMesh = ~doMesh;
            KeyIsDown=0;
        end

        if ( KeyIsDown==1 && KeyCode(quitkey)==1 )
            break;
        end
    else
        oldKeyIsDown = KeyIsDown;
    end

    % Update frame animation counter:
    framecount = framecount + 1;

    % We're done for this frame:

    % Show rendered image 'waitframes' refreshes after the last time
    % the display was updated and in sync with vertical retrace:
    vbl = Screen('Flip', win, vbl + (waitframes - 0.5) * ifi);
    %Screen('Flip', win, 0, 0, 2);
end

% End of display loop: Shutdown actions...
vbl = Screen('Flip', win);

% Calculate and display average framerate:
fps = framecount / (vbl - tstart) %#ok<NOPRT,NASGU>

% Delete last mesh:
if ~isempty(mymesh)
    PsychKinect('DeleteObject', win, mymesh);
    mymesh = [];
end

% Stop kinects capture operation:
PsychKinect('Stop', kinect);

% Close device, release it into standby mode:
PsychKinect('Close', kinect);

% Shutdown high-level driver:
PsychKinect('Shutdown');

% Close onscreen window and release all other ressources:
Screen('CloseAll');

% Reenable Synctests after this simple demo:
Screen('Preference','SkipSyncTests',1);

% Well done!
return

% drawShape does the actual drawing:
function drawShape(mymesh, doMesh, ang, tilt)
% GL needs to be defined as "global" in each subfunction that
% executes OpenGL commands:
global GL
global win

% Backup modelview matrix:
glPushMatrix;

% Setup rotation around axis:
glTranslatef(0,0,-10);
glRotated(ang,0,1,0);
glRotated(tilt,1,0,0);

% Some other rotations/translations to make it look nice:
glTranslatef(0,0,+10);
glRotated(180,0,1,0);
glRotated(180,0,0,1);

% Scale object by a factor of a:
a=10;
glScalef(a,a,a);

% Render!
PsychKinect('RenderObject', win, mymesh, doMesh);

% Restore modelview matrix:
glPopMatrix;

% Done, return to main:
return;

