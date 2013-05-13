function ContrastModulatedNoiseTheElegantStyleDemo(noisesize, staticnoise)
% ContrastModulatedNoiseTheElegantStyleDemo([noisesize=512] [, staticnoise=0])
%
% This demo shows how to render contrast modulated noise efficiently on
% current state of the art graphics hardware by use of the Psychtoolbox
% imaging pipeline and high precision framebuffers.
%
% A rectangular image with random noise of size
% 'noisesize' by 'noisesize' pixels is created, then drawn to the display.
%
% We use alpha-blending tricks to modulate the contrast of the noise in
% realtime. In this demo you can move the mouse pointer to drag around a
% "modulation disk" of 50 pixels diameter. Noise inside the disk has a
% different contrast 'fgcontrast' from the rest of the image, which has a
% contrast of 'bgcontrast' (see top of code for fgcontrast and bgcontrast).
%
% You can press the cursor left and cursor right keys to decrement or
% increment the contrast inside the disk in steps of 0.01 units. Press
% cursor up key to set the inside disk contrast equal to outside disk
% contrast. Press ESCAPE key to finish the demo.
%
% The optional 'staticnoise' flag - if set to 1 - will render a static
% noise image instead of one that changes at each frame. That's faster
% because one doesn't need to recreate the noise texture each frame.
%
% How this works? Basically we use standard Screen 2D drawing commands to
% draw a "contrast values weight map" into the alpha-channel, so the
% alpha-channel encodes contrast values between 0.0 and 1.0.
% Then we draw a noise texture of fixed contrast into the framebuffer
% and the alpha-blending hardware takes care of modulating the
% contrast of the drawn noise texture with the values from the
% alpha-channel.
%
% Hardware requirements: ATI Radeon X-1000 or later, NVidia Geforce 6000 or
% later. Recommended is Radeon HD2000/3000/... or Geforce-8000/9000/...
% hardware for maximum fun!
%
% The allowable contrast values for standard 2D drawing commands are
% limited to the range 0.0 to 1.0 on MacOS/X, they are not limited on
% MS-Windows or GNU/Linux. Precision is set to 16bpc float, ie. about 3
% digits behind the decimal point or 1024 discriminable levels. On the most
% recent ATI HD-2000 and NVidia Geforce 8000 hardware one can lift this limit
% to 32bpc float -- 6.5 digits or about 8 million discriminable levels by a
% one-line code change in this script ;-)

% History:
% 5.9.2007 Written (MK).

% Running on PTB-3?
AssertOpenGL;

if nargin < 1
    noisesize = [];
end

if isempty(noisesize)
    noisesize = 512;
end

if nargin < 2
    staticnoise = [];
end

if isempty(staticnoise)
    staticnoise = 0;
end

bgcontrast = 0.2;
fgcontrast = 0.4;

% Assign control keys: Cursor Left and Right for increasing/decreasing
% contrast of the noise:
KbName('UnifyKeyNames');
leftArrow = KbName('LeftArrow');
rightArrow = KbName('RightArrow');
upArrow = KbName('UpArrow');
escape = KbName('ESCAPE');

% Open onscreen window on screen with maximum id:
screenid=max(Screen('Screens'));

try
    % Open onscreen window: We request a 32 bit per color component
    % floating point framebuffer if it supports alpha-blendig. Otherwise
    % the system shall fall back to a 16 bit per color component
    % framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    [win, winRect] = PsychImaging('OpenWindow', screenid);

    % We use a normalized color range from now on. All color values are
    % specified as numbers between 0.0 and 1.0, instead of the usual 0 to
    % 255 range. This is more intuitive:
    Screen('ColorRange', win, 1, 0);
    
    % Fill the whole onscreen window with a neutral 50% intensity
    % background color and an alpha channel value of 'bgcontrast'.
    % This becomes the clear color. After each Screen('Flip'), the
    % backbuffer will be cleared to this neutral 50% intensity gray
    % and a default 'bgcontrast' background noise contrast level:
    Screen('FillRect', win, [0.5 0.5 0.5 bgcontrast]);
    
    i=0;
    tonset = [];

    % Create first noisematrix outside loop:
    noisematrix = randn(noisesize);
    
    % Initially place the mouse cursor at the center of the window:
    [x, y] = RectCenter(winRect);
    SetMouse(x, y, win);
    HideCursor;
    
    % Initially sync us to retrace:
    Screen('Flip', win);
    
    % Main stimulus drawing loop:
    while 1
        % Increment framecounter:
        i=i+1;
        
        % If staticnoise is set to 1 then we only generate the noise
        % texture once -- This way we can disentangle "noise creation
        % overhead" from actual drawing overhead:
        if i==1 || staticnoise == 0
            % Convert Matlab 'noisematrix' to 16 bpc floating point noise texture:
            noisetex = Screen('MakeTexture', win, noisematrix, [], [], 1);
        end
        
        % Disable alpha-blending, so we can just overwrite the framebuffer
        % with our new pixels:
        Screen('Blendfunction', win, GL_ONE, GL_ZERO);
        
        % Now we overdraw some regions of the onscreen windows alpha-channel
        % with our "modulation" image - a image that contains alpha values
        % which encode a different contrast 'fgcontrast'. After this drawing op,
        % the alpha-channel will contain the final "contrast modulation landscape":
        Screen('DrawDots', win, [x y], 50, [0.5 0.5 0.5 fgcontrast], [], 1);

        % On some graphics hardware + operating system combos, 'DrawDots'
        % is not able to draw round dots when Screen('ColorRange', win, 1,
        % 1); is used. In such cases, just use 'FillOval' to draw round
        % dots:
        % Screen('FillOval', win, [0.5 0.5 0.5 fgcontrast], OffsetRect([0 0 50 50], x-25, y-25));

        % Now we draw the noise texture and use alpha-blending of
        % the drawn noise color pixels with the destination alpha-channel,
        % thereby multiplying the incoming color values with the stored
        % alpha values -- effectively a contrast modulation. The GL_ONE
        % means that we add the final contrast modulated noise pixels to
        % the current content of the window == the neutral gray background.
        Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE);
        
        % The extra zero at the end forcefully disables bilinear filtering. This is
        % not strictly neccessary on correctly working hardware, but an extra
        % precaution to make sure that the noise values are blitted
        % one-to-one into the offscreen window:
        Screen('DrawTexture', win, noisetex, [], [], [], 0);
            
        % At this point, the final image should be ready in the backbuffer
        % of our onscreen window. Ready to flip it onscreen...
        
        % Cleanup: Release our noise texture in the dynamic-noise case were
        % we recreate one for each drawn frame:
        if ~staticnoise
            Screen('Close', noisetex);
        end

        % Tell PTB that all drawing commands are done now. This allows
        % the graphics hardware to perform all drawing and image processing
        % in parallel while we execute Matlab code for non-graphics related
        % stuff, in our case random noise creation, keyboard and mouse
        % queries:
        Screen('DrawingFinished', win);
                
        % Now all the non-Screen() stuff:
        if ~staticnoise
            % Generate a noisesize x noisesize matrix of random noise with mean zero,
            % stddev 1.0 for use in the next loop iteration:
            noisematrix = randn(noisesize);
        end
        
        % Keyboard queries:
        [isdown secs keycode] = KbCheck;
        if isdown
            if keycode(escape)
                break;
            end
            
            if keycode(leftArrow)
                fgcontrast = max(0, fgcontrast - 0.01);
            end
            
            if keycode(rightArrow)
                fgcontrast = min(1, fgcontrast + 0.01);
            end
            
            if keycode(upArrow)
                fgcontrast = bgcontrast;
            end
        end

        % Query mouse position: This will be the center of our "modulation disk":
        [x, y] = GetMouse(win);

        % Ready. Request stimulus onset:
        tonset(i) = Screen('Flip', win);

        % Ready. Next loop iteration:
    end

    % One final flip:
    Screen('Flip', win);

    % Done. Close screen and finish:
    ShowCursor;
    Screen('CloseAll');
    
    % Compute avg. computation time for redraw:
    avgredrawtime = mean(diff(tonset)) * 1000
    %plot(diff(tonset));
    
    % Done.
    return;
catch
    % Error. Close screen, show cursor, rethrow error:
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
