function StereoDemo(stereoMode)
% StereoDemo(stereoMode)
%
% Demo on how to use OpenGL-Psychtoolbox to present stereoscopic stimuli.
%
% DEPRECATED! THIS DEMO ONLY DEMONSTRATES THE LEGACY STEREO SUPPORT!
% If you have a recent graphics card, have a look at ImagingStereoDemo.m
% instead. It demonstrates a much better, more flexible, robust way of
% presenting stereo stimuli with Psychtoolbox.
%
%
% Press any key to abort demo any time.
%
% stereoMode specifies the type of stereo display algorithm to use:
%
% 0 == Mono display - No stereo at all.
%
% 1 == Flip frame stereo (temporally interleaved) - You'll need shutter
% glasses that are supported by the operating system, e.g., the
% CrystalEyes-Shutterglasses.
%
% 2 == Top/bottom image stereo with lefteye=top also for use with special
% CrystalEyes-hardware.
%
% 3 == Same, but with lefteye=bottom.
%
% 4 == Free fusion (lefteye=left, righteye=right): This - together wit a
% screenid of zero - is what you'll want to use on MS-Windows dual-display
% setups for stereo output.
%
% 5 == Cross fusion (lefteye=right ...)
%
% 6-9 == Different modes of anaglyph stereo for color filter glasses:
%
% 6 == Red-Green
% 7 == Green-Red
% 8 == Red-Blue
% 9 == Blue-Red
%
% 10 == Dual-Window stereo: Open two onscreen windows, first one will
% display left-eye view, 2nd one right-eye view. Direct all drawing and
% flip commands to the first window, PTB will take care of the rest. This
% mode is mostly useful for dual-display stereo on MacOS/X. It only works
% on reasonably modern graphics hardware, will abort with an error on
% unsupported hardware.
%
% 11 == Like mode 1 (frame-sequential) but using Screen's built-in method,
% instead of the native method supported by your graphics card.
%

% Authors:
% Finnegan Calabro  - fcalabro@bu.edu
% Mario Kleiner     - mario.kleiner at tuebingen.mpg.de
%

if nargin < 1
    % When no stereoMode is given, we start with Anaglyph stereo, a safe
    % choice on any display hardware...
    stereoMode=8;
end;

% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL;

try
    % Get the list of Screens and choose the one with the highest screen number.
    % Screen 0 is, by definition, the display with the menu bar. Often when
    % two monitors are connected the one without the menu bar is used as
    % the stimulus display.  Chosing the display with the highest dislay number is
    % a best guess about where you want the stimulus displayed.
    scrnNum = max(Screen('Screens'));

    % Windows-Hack: If mode 4 or 5 is requested, we select screen zero
    % as target screen: This will open a window that spans multiple
    % monitors on multi-display setups, which is usually what one wants
    % for this mode.
    if IsWin && (stereoMode==4 || stereoMode==5)
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
    
    % Stimulus settings:
    numDots = 1000;
    vel = 1;   % pix/frames
    dotSize = 4;
    dots = zeros(3, numDots);

    xmax = 300;
    ymax = xmax;

    f = 4*pi/xmax;
    amp = 16;

    dots(1, :) = 2*(xmax)*rand(1, numDots) - xmax;
    dots(2, :) = 2*(ymax)*rand(1, numDots) - ymax;

    % Open double-buffered onscreen window with the requested stereo mode:
    [windowPtr, windowRect]=Screen('OpenWindow', scrnNum, BlackIndex(scrnNum), [], [], [], stereoMode);
    
    % This call would change the color of blue-sync lines in
    % frame-sequential stereo mode to blue:
    % SetStereoBlueLineSyncParameters(windowPtr, [], [], [0,0,1]);
    
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

    % Perform a flip to sync us to vbl and take start-timestamp in t:
    t = Screen('Flip', windowPtr);
    onset = t;
    
    % Run until a key is pressed:
    while ~KbCheck
        % Compute dot positions and offsets for this frame:
        center = center + [xvel yvel];
        if center(1) > xmax || center(1) < -xmax
            xvel = -xvel;
        end

        if center(2) > ymax || center(2) < -ymax
            yvel = -yvel;
        end

        dots(3, :) = -amp.*exp(-(dots(1, :) - center(1)).^2 / (2*sigma*sigma)).*exp(-(dots(2, :) - center(2)).^2 / (2*sigma*sigma));

        % Select left-eye image buffer for drawing:
        Screen('SelectStereoDrawBuffer', windowPtr, 0);
        % Draw left stim:
        Screen('DrawDots', windowPtr, dots(1:2, :) + [dots(3, :)/2; zeros(1, numDots)], dotSize, col1, [windowRect(3:4)/2], 1);

        % Select right-eye image buffer for drawing:
        Screen('SelectStereoDrawBuffer', windowPtr, 1);
        % Draw right stim:
        Screen('DrawDots', windowPtr, dots(1:2, :) - [dots(3, :)/2; zeros(1, numDots)], dotSize, col2, [windowRect(3:4)/2], 1);

        % Take timestamp of stimulus-onset after displaying the new stimulus
        % and record it in vector t:
        onset = Screen('Flip', windowPtr, onset+0.001);
        t = [t onset];
    end

    % Done. Close the onscreen window:
    Screen('CloseAll')

    % Compute and show timing statistics:
    dt = t(2:end) - t(1:end-1);
    disp(sprintf('N.Dots\tMean (s)\tMax (s)\t%%>20ms\t%%>30ms\n'));
    disp(sprintf('%d\t%5.3f\t%5.3f\t%5.2f\t%5.2f\n', numDots, mean(dt), max(dt), sum(dt > 0.020)/length(dt), sum(dt > 0.030)/length(dt)));

    % We're done.

catch
    % Executes in case of an error: Closes onscreen window:
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end;
