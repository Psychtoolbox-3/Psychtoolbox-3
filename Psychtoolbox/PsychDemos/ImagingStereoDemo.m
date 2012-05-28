function ImagingStereoDemo(stereoMode, usedatapixx, writeMovie)
% ImagingStereoDemo([stereoMode=8][, usedatapixx = 0][, writeMovie = 0])
%
% Demo on how to use OpenGL-Psychtoolbox to present stereoscopic stimuli
% when the Psychtoolbox imaging pipeline is enabled. Use of the imaging
% pipeline allows for more flexible and high quality stereo display modes,
% but it requires graphics hardware with support for at least framebuffer
% objects and Shadermodel 2.0. See the Psychtoolbox Wiki about gfx-hardware
% recommendations. The demo also shows how to configure the pipeline to
% restrict image processing to some subregion of the display, e.g., to save
% some computation time on low-end hardware.
%
% Press escape key to abort demo, space key to toggle modes of specific
% algorithms.
%
% 'stereoMode' specifies the type of stereo display algorithm to use:
%
% 0 == Mono display - No stereo at all.
%
% 1 == Flip frame stereo (temporally interleaved) - You'll need shutter
% glasses that are supported by the operating system, e.g., the
% CrystalEyes-Shutterglasses. Psychtoolbox will automatically generate blue
% sync lines at the bottom of the display.
%
% 2 == Top/bottom image stereo with lefteye=top also for use with special
% CrystalEyes-hardware.
%
% 3 == Same, but with lefteye=bottom.
%
% 4 == Free fusion (lefteye=left, righteye=right): Left-eye view displayed
% in left half of window, right-eye view displayed in right-half of window.
% Use this for dual-display setups (binocular video goggles, haploscopes,
% polarized stereo setups etc.)
%
% 5 == Cross fusion (lefteye=right ...): Like mode 4, but with views
% switched.
%
% 6-9 == Different modes of anaglyph stereo for color filter glasses:
%
% 6 == Red-Green
% 7 == Green-Red
% 8 == Red-Blue
% 9 == Blue-Red
%
% 10 == Like mode 4, but for use on Mac OS/X with dualhead display setups.
%
% 11 == Like mode 1 (frame-sequential) but using Screen's built-in method,
% instead of the native method supported by your graphics card.
%
% 100 == Interleaved line stereo: Left eye image is displayed in even
% scanlines, right eye image is displayed in odd scanlines.
%
% 101 == Interleaved column stereo: Left eye image is displayed in even
% columns, right eye image is displayed in odd columns. Typically used for
% auto-stereoscopic displays, e.g., lenticular sheet or parallax barrier
% displays.
%
% If you have a different set of filter glasses, e.g., red-magenta, you can
% simply select one of above modes, then use the
% SetStereoAnglyphParameters() command to change color gain settings,
% thereby implementing other anaglyph color combinations.
%
% 'usedatapixx' If provided and set to a non-zero value, will setup a
% connected VPixx DataPixx device for stereo display.
%
% 'writeMovie' If provided and set to a non-zero value will write a
% Quicktime movie file 'MyTestMovie.mov' into the current working directory
% which captures the full performance of this demo. A setting of 1 will
% only write video, a setting of 2 will also write an audio track with
% a sequence of ten successive beep tones of 1 sec duration.
%
% Authors:
% Finnegan Calabro  - fcalabro@bu.edu
% Mario Kleiner     - mario.kleiner at tuebingen.mpg.de
%

% We start of with non-inverted display:
inverted = 0;

% Default to stereoMode 8 -- Red-Green stereo:
if nargin < 1
    stereoMode = [];
end

if isempty(stereoMode)
    stereoMode = 8;
end;

if nargin < 2
    usedatapixx = [];
end

if isempty(usedatapixx)
    usedatapixx = 0;
end

if nargin < 3
    writeMovie = [];
end

if isempty(writeMovie)
    writeMovie = 0;
end

% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL;

% Define response key mappings, unify the names of keys across operating
% systems:
KbName('UnifyKeyNames');
space = KbName('space');
escape = KbName('ESCAPE');

%try
% Get the list of Screens and choose the one with the highest screen number.
% Screen 0 is, by definition, the display with the menu bar. Often when
% two monitors are connected the one without the menu bar is used as
% the stimulus display.  Chosing the display with the highest dislay number is
% a best guess about where you want the stimulus displayed.
scrnNum = max(Screen('Screens'));

% Increase level of verbosity for debug purposes:
%Screen('Preference', 'Verbosity', 6);

% Windows-Hack: If mode 4 or 5 is requested, we select screen zero
% as target screen: This will open a window that spans multiple
% monitors on multi-display setups, which is usually what one wants
% for this mode.
if IsWin & (stereoMode==4 | stereoMode==5)
    scrnNum = 0;
end

% Dual display dual-window stereo requested?
if stereoMode == 10
    % Yes. Do we have at least two separate displays for both views?
    if length(Screen('Screens')) < 2
        error('Sorry, for stereoMode 10 you''ll need at least 2 separate display screens in non-mirrored mode.');
    end

    if ~IsWin
        % Assign left-eye view (the master window) to main display:
        scrnNum = 0;
    else
        % Assign left-eye view (the master window) to main display:
        scrnNum = 1;
    end
end

% Open double-buffered onscreen window with the requested stereo mode,
% setup imaging pipeline for additional on-the-fly processing:

% Prepare pipeline for configuration. This marks the start of a list of
% requirements/tasks to be met/executed in the pipeline:
PsychImaging('PrepareConfiguration');

if usedatapixx
    % Tell PTB we want to display on a DataPixx device:
    PsychImaging('AddTask', 'General', 'UseDataPixx');
end

% Ask to restrict stimulus processing to some subarea (ROI) of the
% display. This will only generate the stimulus in the selected ROI and
% display the background color in all remaining areas, thereby saving
% some computation time for pixel processing: We select the center
% 512x512 pixel area of the screen:
if ~ismember(stereoMode, [100, 101])
    PsychImaging('AddTask', 'AllViews', 'RestrictProcessing', CenterRect([0 0 512 512], Screen('Rect', scrnNum)));
end

% stereoMode 100 triggers scanline interleaved display:
if stereoMode == 100
    PsychImaging('AddTask', 'General', 'InterleavedLineStereo', 0);
end

% stereoMode 101 triggers column interleaved display:
if stereoMode == 101
    PsychImaging('AddTask', 'General', 'InterleavedColumnStereo', 0);
end

% Consolidate the list of requirements (error checking etc.), open a
% suitable onscreen window and configure the imaging pipeline for that
% window according to our specs. The syntax is the same as for
% Screen('OpenWindow'):
[windowPtr, windowRect]=PsychImaging('OpenWindow', scrnNum, 0, [], [], [], stereoMode);

% Oldstyle:   [windowPtr, windowRect]=Screen('OpenWindow', scrnNum, 0, [], [], [], stereoMode);

if stereoMode == 10
    % In dual-window, dual-display mode, we open the slave window on
    % the secondary screen. Please note that, after opening this window
    % with the same parameters as the "master-window", we won't touch
    % it anymore until the end of the experiment. PTB will take care of
    % managing this window automatically as appropriate for a stereo
    % display setup. That is why we are not even interested in the window
    % handles of this window:
    if IsWin
        slaveScreen = 2;
    else
        slaveScreen = 1;
    end
    Screen('OpenWindow', slaveScreen, BlackIndex(slaveScreen), [], [], [], stereoMode);
end

% Stimulus settings:
numDots = 1000;
vel = 1;   % pix/frames
dotSize = 8;
dots = zeros(3, numDots);

xmax = RectWidth(windowRect)/2;
ymax = RectHeight(windowRect)/2;
if stereoMode == 100
    xmax = xmax/4;
    ymax = ymax/2;
else
    xmax = min(xmax, ymax) / 2;
    ymax = xmax;
end

f = 4*pi/xmax;
amp = 16;

dots(1, :) = 2*(xmax)*rand(1, numDots) - xmax;
dots(2, :) = 2*(ymax)*rand(1, numDots) - ymax;

% Set color gains. This depends on the anaglyph mode selected. The
% values set here as default need to be fine-tuned for any specific
% combination of display device, color filter glasses and (probably)
% lighting conditions and subject. The current settings do ok on a
% MacBookPro flat panel.
switch stereoMode
    case 6,
        SetAnaglyphStereoParameters('LeftGains', windowPtr,  [1.0 0.0 0.0]);
        SetAnaglyphStereoParameters('RightGains', windowPtr, [0.0 0.6 0.0]);
    case 7,
        SetAnaglyphStereoParameters('LeftGains', windowPtr,  [0.0 0.6 0.0]);
        SetAnaglyphStereoParameters('RightGains', windowPtr, [1.0 0.0 0.0]);
    case 8,
        SetAnaglyphStereoParameters('LeftGains', windowPtr, [0.4 0.0 0.0]);
        SetAnaglyphStereoParameters('RightGains', windowPtr, [0.0 0.2 0.7]);
    case 9,
        SetAnaglyphStereoParameters('LeftGains', windowPtr, [0.0 0.2 0.7]);
        SetAnaglyphStereoParameters('RightGains', windowPtr, [0.4 0.0 0.0]);
    otherwise
        %error('Unknown stereoMode specified.');
end

% Initially fill left- and right-eye image buffer with black background
% color:
Screen('SelectStereoDrawBuffer', windowPtr, 0);
Screen('FillRect', windowPtr, BlackIndex(scrnNum));
Screen('SelectStereoDrawBuffer', windowPtr, 1);
Screen('FillRect', windowPtr, BlackIndex(scrnNum));

% Show cleared start screen:
Screen('Flip', windowPtr);

% Set up alpha-blending for smooth (anti-aliased) drawing of dots:
Screen('BlendFunction', windowPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

col1 = WhiteIndex(scrnNum);
col2 = col1;
i = 1;
keyIsDown = 0;
center = [0 0];
sigma = 50;
xvel = 2*vel*rand(1,1)-vel;
yvel = 2*vel*rand(1,1)-vel;

Screen('Flip', windowPtr);

% Maximum number of animation frames to show:
nmax = 100000;

% Optionally create a Quicktime movie file 'MyTestMovie.mov' in the
% current directory. The file will record a movie of this performance
% with video frames of size 512 x 512 pixels at a framerate of 60fps.
if writeMovie
    if writeMovie > 1
        if ~IsOSX(1)
            % Add a sound track to the movie: 2 channel stereo, 48 kHz:
            movie = Screen('CreateMovie', windowPtr, ['MyTestMovie.mov'], 512, 512, 30, ':CodecSettings=AddAudioTrack=2@48000');
        else
            % Same, but for 64-Bit OSX, where the default codec does not
            % work due to too old GStreamer version. Use XVid instead:
            movie = Screen('CreateMovie', windowPtr, ['MyTestMovie.avi'], 512, 512, 30, ':CodecSettings=AddAudioTrack=2@48000:CodecType=VideoCodec=xvidenc profile=244 max-key-interval=10 bitrate=9708400 quant-type=1');
        end
        % Create a sequence of 10 tones, each of 1 second duration, each 100 Hz higher
        % than its predecessor. Each of the two stereo channels gets a slightly different sound:
        for freq=100:100:1000
            Screen('AddAudioBufferToMovie', movie, [0.8 * MakeBeep(freq, 1, 48000); 0.8 * MakeBeep(freq*1.2, 1, 48000)]);
        end
        nmax = 300;
    else
        % Only video, no sound:
        movie = Screen('CreateMovie', windowPtr, ['MyTestMovie.mov'], 512, 512, 30);
    end

    % Other examples of codec settings:
    %
    %        movie = Screen('CreateMovie', windowPtr, 'WinXPTest.avi', 640, 480, 30, ':CodecType=VideoCodec=x264enc speed-preset=5 key-int-max=30 bitrate=20000 profile=3');
    %        movie = Screen('CreateMovie', windowPtr, 'WinXPTest.avi', 640, 480, 30, ':CodecType=VideoCodec=x264enc speed-preset=5 bitrate=20000 profile=3');
    %        movie = Screen('CreateMovie', windowPtr, 'WinXPTest.avi', 640, 480, 30, ':CodecType=theoraenc');
    %        movie = Screen('CreateMovie', windowPtr, 'WinXPTest.avi', 640, 480, 30, ':CodecType=theoraenc AddAudioTrack');
    %        movie = Screen('CreateMovie', windowPtr, 'WinXPTest.avi', 320, 240, 30, ':CodecType=VideoCodec=xvidenc profile=244 max-key-interval=10 bitrate=9708400 quant-type=1');
    %        A full gst-launch style gstreamer launch line with everything
    %        manually controlled. This for Flash video encoding (.flv):
    %        movie = Screen('CreateMovie', windowPtr, ['MyTestMovie.flv'], 512, 512, 30, 'gst-launch appsrc name=ptbvideoappsrc do-timestamp=0 stream-type=0 max-bytes=0 block=1 is-live=0 emit-signals=0 ! capsfilter caps="video/x-raw-rgb, bpp=(int)32, depth=(int)32, endianess=(int)4321, alpha_mask=(int)-16777216, red_mask=(int)16711680, green_mask=(int)65280, blue_mask=(int)255, width=(int)512, height=(int)512, framerate=30/1" ! videorate ! ffmpegcolorspace ! ffenc_flv ! ffmux_flv ! filesink name=ptbfilesink async=0 location="MyTestMovie.flv"');
end

% Perform a flip to sync us to vbl and take start-timestamp in t:
t = Screen('Flip', windowPtr);

% Run until a key is pressed:
while length(t) < nmax

    % Select left-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', windowPtr, 0);

    % Draw left stim:
    Screen('DrawDots', windowPtr, dots(1:2, :) + [dots(3, :)/2; zeros(1, numDots)], dotSize, col1, [windowRect(3:4)/2], 1);
    Screen('FrameRect', windowPtr, [255 0 0], [], 5);

    % Select right-eye image buffer for drawing:
    Screen('SelectStereoDrawBuffer', windowPtr, 1);

    % Draw right stim:
    Screen('DrawDots', windowPtr, dots(1:2, :) - [dots(3, :)/2; zeros(1, numDots)], dotSize, col2, [windowRect(3:4)/2], 1);
    Screen('FrameRect', windowPtr, [0 255 0], [], 5);

    % Tell PTB drawing is finished for this frame:
    Screen('DrawingFinished', windowPtr);

    % Now all non-drawing tasks:

    % Compute dot positions and offsets for next frame:
    center = center + [xvel yvel];
    if center(1) > xmax | center(1) < -xmax
        xvel = -xvel;
    end

    if center(2) > ymax | center(2) < -ymax
        yvel = -yvel;
    end

    dots(3, :) = -amp.*exp(-(dots(1, :) - center(1)).^2 / (2*sigma*sigma)).*exp(-(dots(2, :) - center(2)).^2 / (2*sigma*sigma));

    % Keyboard queries and key handling:
    [pressed dummy keycode] = KbCheck;
    if pressed
        % SPACE key toggles between non-inverted and inverted display:
        if keycode(space) & ismember(stereoMode, [6 7 8 9]);
            while KbCheck; end;
            inverted = 1 - inverted;
            if inverted
                % Set inverted mode:
                SetAnaglyphStereoParameters('InvertedMode', windowPtr);
            else
                % Set standard mode:
                SetAnaglyphStereoParameters('StandardMode', windowPtr);
            end
        end

        % ESCape key exits the demo:
        if keycode(escape)
            break;
        end
    end

    % Add a screenshot of the center 512 x 512 pixels as a new video frame to the movie file, if any:
    if writeMovie
        Screen('AddFrameToMovie', windowPtr, CenterRect([0 0 512 512], Screen('Rect', scrnNum)), 'backBuffer');
    end

    % Flip stim to display and take timestamp of stimulus-onset after
    % displaying the new stimulus and record it in vector t:
    onset = Screen('Flip', windowPtr);
    t = [t onset];
end

% Last Flip:
Screen('Flip', windowPtr);

% Finalize and close movie file, if any:
if writeMovie
    Screen('FinalizeMovie', movie);
end

% Done. Close the onscreen window:
Screen('CloseAll')

% Compute and show timing statistics:
dt = t(2:end) - t(1:end-1);
disp(sprintf('N.Dots\tMean (s)\tMax (s)\t%%>20ms\t%%>30ms\n'));
disp(sprintf('%d\t%5.3f\t%5.3f\t%5.2f\t%5.2f\n', numDots, mean(dt), max(dt), sum(dt > 0.020)/length(dt), sum(dt > 0.030)/length(dt)));

% We're done.
return;
%catch
% Executes in case of an error: Closes onscreen window:
Screen('CloseAll');
psychrethrow(psychlasterror);
%end;
