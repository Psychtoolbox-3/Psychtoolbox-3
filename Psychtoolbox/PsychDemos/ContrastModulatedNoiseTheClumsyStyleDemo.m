function ContrastModulatedNoiseTheClumsyStyleDemo(noisesize, staticnoise)
% ContrastModulatedNoiseTheClumsyStyleDemo([noisesize=512] [, staticnoise=0])
%
% This demo shows how to render contrast modulated noise with old, outdated
% graphics hardware. A rectangular image with random noise of size
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
% alpha-channel encodes contrast values between 0.0 and 1.0 in steps of
% 1/256th. Then we draw a noise texture of fixed contrast into the
% framebuffer and the alpha-blending hardware takes care of modulating the
% contrast of the drawn noise texture with the values from the
% alpha-channel.
%
% Hardware requirements: ATI Radeon 8500 or later, NVidia Geforce2 or later
% on MacOS/X. Not Intel-GMA 950 of Intel MacBook!
%
% Limitations: This demo is meant to run on outdated graphics hardware, so
% it can't make use of the Psychtoolbox imaging pipeline, which would allow
% for a very elegant, flexible and fast implementation. We use a multi-pass
% approach here to work around hardware restrictions as much as possible.
% Some of the limitations of this approach: Contrast modulation is
% restricted to 0.0 - 1.0 and accuracy is limited to 256 levels. Clamping
% of your noise matrix (has to fit in the "normal" 8 bit 256 gray level
% range) can change the noise distribution in unwanted ways. The multi-pass
% approach used here is pretty inefficient, so redraw rates are pretty low,
% e.g., max. 15 fps on a fast MacBookPro with ATI Radeon X1600.
%
% If you have recent graphics hardware (ATI Radeon X1000 or later, NVidia
% Geforce 6000 and later), take a look at ContrastModulatedNoiseTheElegantStyleDemo.m
% instead. That demo takes advantage of the imaging pipeline, greatly
% simplifying stimulus programming and lifting basically all restrictions
% of the approach presented here. Achievable redraw rates are also much
% higher that way!

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

% Open onscreen window on screen with maximum id: The window gets a
% background clear color of RGB==128==neutral 50% grey and a alpha-value of
% zero. The alpha channel will be used to store the ...
screenid=max(Screen('Screens'));

try

    % Initialize OpenGL, but only for 2D drawing. We need this to be able
    % to use the low-level OpenGL command glBlendFunc():
    InitializeMatlabOpenGL([],[],1);
    [win, winRect] = Screen('OpenWindow', screenid, 128);

    % Create an offscreen window for noise with values >= 0:
    positivenoise = Screen('OpenOffscreenWindow', win, [0 0 0 0]);

    % Create an offscreen window for noise with values < 0:
    negativenoise = Screen('OpenOffscreenWindow', win, [0 0 0 0]);

    i=0;
    tonset = [];

    % Initially place the mouse cursor at the center of the window:
    [x, y] = RectCenter(winRect);
    SetMouse(x, y, win);
    HideCursor;
    
    % Main stimulus drawing loop:
    while 1
        % Increment framecounter:
        i=i+1;
        
        % Query mouse position: This will be the center of our "modulation disk":
        [x, y] = GetMouse(win);
        
        % Compute target rectangle for our drawn modulation oval:
        dstrect = CenterRectOnPoint([0 0 50 50], x, y);

        % If staticnoise is set to 1 then we only generate the noise
        % textures once -- This way we can disentangle "noise creation
        % overhead" from actual drawing overhead:
        if i==1 || staticnoise == 0
            % Generate a noisesize x noisesize matrix of random noise, mean zero,
            % stddev 1.0:
            noisematrix = randn(noisesize) * 255;

            % The following conversion has a problem: We map the magnitudes
            % between zero and 1.0 to 0-255. Values outside that magnitude
            % range will create pixel values > 255 and will be therefore
            % wrapped around (remapped) to the range 0-255. This will alter
            % the noise distribution considerably! Alternatively one could
            % generate noise with a way smaller standard deviation to avoid
            % this in most cases...
            
            % Convert to texture which only contains the positive values:
            positivenoisetex = Screen('MakeTexture', win, max(0, noisematrix));

            % Convert to texture which only contains the negative values, but
            % mirrored to become positive values:
            negativenoisetex = Screen('MakeTexture', win, -1 * min(0, noisematrix));
        end
            
        % Enable standard additive blending:
        glBlendEquation(GL.FUNC_ADD);

        % We use the alpha-channels of the offscreen windows to define the
        % modulationgains for the noise. Example: Alpha value zero = contrast
        % zero. Alpha value 255 = Maximum contrast of 1.0, Alpha value 128 =
        % Contrast of 128/255 == about 50% or 0.5.
        % First we need to initialize the whole windows alpha-channel with the
        % default "background" contrast bgcontrast. We disable alpha-blending
        % for this purpose:

        Screen('Blendfunction', positivenoise, GL_ONE, GL_ZERO);
        Screen('FillRect', positivenoise, [0 0 0 255 * bgcontrast]);

        % Now we overdraw some regions of the alpha-channel with our "modulation"
        % image - a image that contains alpha values which encode a different
        % contrast 'fgcontrast'. After this drawing op, the alpha-channel will
        % contain the final "contrast modulation landscape":
        %Screen('DrawDots', positivenoise, [x y], 50, [0 0 0 255 * fgcontrast], [], 1);
        Screen('FillOval', positivenoise, [0 0 0 255 * fgcontrast], dstrect);

        % Now we draw the positivenoise texture and use alpha-blending of
        % the drawn noise color pixels with the destination alpha-channel,
        % thereby multiplying the incoming color values with the stored
        % alpha values -- effectively a contrast modulation:
        Screen('Blendfunction', positivenoise, GL_DST_ALPHA, GL_ZERO);
        
        % The extra zero at the end forcefully disables bilinear filtering. This is
        % not strictly neccessary on correctly working hardware, but an extra
        % precaution to make sure that the noise values are blitted
        % one-to-one into the offscreen window:
        Screen('DrawTexture', positivenoise, positivenoisetex, [], [], [], 0);
        
        % Now the positivenoise offscreen window contains a correctly
        % contrast modulated noise image - but only the positive noise
        % values. We now repeat the procedure for the negative noise
        % values:
        
        Screen('Blendfunction', negativenoise, GL_ONE, GL_ZERO);
        Screen('FillRect', negativenoise, [0 0 0 255 * bgcontrast]);

        % Now we overdraw some regions of the alpha-channel with our "modulation"
        % image - a image that contains alpha values which encode a different
        % contrast 'fgcontrast'. After this drawing op, the alpha-channel will
        % contain the final "contrast modulation landscape":
        % Screen('DrawDots', negativenoise, [x y], 50, [0 0 0 255 * fgcontrast], [], 1);
        Screen('FillOval', negativenoise, [0 0 0 255 * fgcontrast], dstrect);

        % Now we draw the positivenoise texture and use alpha-blending of
        % the drawn noise color pixels with the destination alpha-channel,
        % thereby multiplying the incoming color values with the stored
        % alpha values -- effectively a contrast modulation:
        Screen('Blendfunction', negativenoise, GL_DST_ALPHA, GL_ZERO);
        
        % The extra zero at the end forcefully disables bilinear filtering. This is
        % not strictly neccessary on correctly working hardware, but an extra
        % precaution to make sure that the noise values are blitted
        % one-to-one into the offscreen window:
        Screen('DrawTexture', negativenoise, negativenoisetex, [], [], [], 0);
        
        % Now for the compositing of the final image into our onscreen
        % windows backbuffer: First draw positive noise image in additive
        % mode -- adding contrast modulated positive noise pixels to the
        % grey 50% background color:
        Screen('Blendfunction', win, GL_ONE, GL_ONE);
        Screen('DrawTexture', win, positivenoise, [], [], [], 0);
        
        % Now draw negative noise image in subtractive mode: Subtract noise
        % pixels from the current content of the framebuffer:
        glBlendEquation(GL.FUNC_REVERSE_SUBTRACT);
        Screen('DrawTexture', win, negativenoise, [], [], [], 0);
        
        % Enable standard additive blending:
        glBlendEquation(GL.FUNC_ADD);
        
        % At this point, the final image should be ready in the backbuffer
        % of our onscreen window. Ready to flip it onscreen:
        tonset(i) = Screen('Flip', win);

        % Release our noise textures in the dynamic-noise case.
        if ~staticnoise
            Screen('Close', positivenoisetex);
            Screen('Close', negativenoisetex);
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

        % Ready. Next loop iteration:
    end

    % Done. Close screen and finish:
    ShowCursor;
    Screen('CloseAll');
    
    % Compute avg. computation time for redraw:
    avgredrawtime = mean(diff(tonset)) * 1000
    plot(diff(tonset));
    
    % Done.
    return;
catch
    % Error. Close screen, show cursor, rethrow error:
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
