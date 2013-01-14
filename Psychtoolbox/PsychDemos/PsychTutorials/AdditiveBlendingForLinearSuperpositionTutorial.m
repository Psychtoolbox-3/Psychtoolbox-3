function AdditiveBlendingForLinearSuperpositionTutorial(outputdevice, overlay, colorclut, doGainCorrection)
% AdditiveBlendingForLinearSuperpositionTutorial([outputdevice='None'] [, overlay=1] [, colorclut=0] [, doGainCorrection=0]);
%
% Illustrates use of floating point textures in combination with
% source-weighted additive alpha blending to create linear superpositions
% of image patches, in this case of two superimposed gratings. Due to the
% use of additive alpha-blending, you'll see that the two gratings will
% superimpose onto each other in a mathematically correct way -- pixel-wise
% addition/subtraction of luminance values. The demo uses a 32 bit floating
% point framebuffer on the latest hardware. This allows for an effective 23
% bits of precision in all math done and in the final stimuli - more than
% any display device in existence could resolve. On previous generation
% hardware (older than NVidia Geforce 88000 or ATI Radeon HD2000), alpha
% blending isn't supported in 32 bpc float precision. Therefore the demo
% will select 16 bpc floating point precision, where alpha blending works.
% This way the effective precision is 11 bits, a bit less than what special
% display devices can resolve. However, final gamma correction is done in
% full 23 bits precision, so you can make use of the extra bits for proper
% display linearization.
%
% The demo also demonstrates Psychtoolbox support for high precision
% display output devices, ie., display extenders, graphics cards and
% display modes that allow for luminance display with more than the
% standard 8 bits precision. By default the standard 8 bit framebuffer of
% standard graphics cards is demonstrated.
%
% However, by providing the first optional parameter 'outputdevice', you
% can select amongst all devices supported by PTB:
%
% 'NoneNoGamma' - Same as 'None': Use standard 8 bit framebuffer, but
% disable the gamma correction provided by PTB's imaging pipeline. This is
% usually not what you want, but it allows to test how much faster the
% display runs without gamma correction.
%
% 'PseudoGray' - PseudoGray display, also known as "Bit stealing". This
% technique allows to create the perception of up to 1786 different
% luminance levels on standard 8 bit graphics hardware by use of some
% clever color rendering trick. See "help CreatePseudoGrayLUT" for
% references and details.
%
% 'Native10Bit' - Enables the native 10 bpc framebuffer support on ATI
% Radeon X1xxx / HDxxx GPU's when used under Linux or OS/X with the
% PsychtoolboxKernelDriver loaded (see "help PsychtoolboxKernelDriver" on
% how to do that). These GPU's do support 10 bits per color channel when
% this special mode is used. If you try this option on MS-Windows, or
% without the driver loaded or with a different GPU, it will just fail.
%
% 'VideoSwitcher' - Enable the Xiangrui Li et al. VideoSwitcher, a special
% type of video attenuator (see "help PsychVideoSwitcher") in standard
% "simple" mode.
%
% 'VideoSwitcherCalibrated' - Enable the Xiangrui Li et al. VideoSwitcher,
% but use the more complex (and more accurate?) mode with calibrated lookup
% tables (see "help PsychVideoSwitcher").
%
% 'Attenuator' - Enable support for standard Pelli & Zhang style video
% attenuators by use of lookup tables.
%
% Then we have support for the different modes of operation of the
% Cambridge Research Systems Bits++ box:
%
% 'Mono++' - Use 14 bit mono output mode, either with color index overlay
% (if the optional 2nd 'overlay' flag is set to 1, which is the default),
% or without color index overlay.
%
% 'Color++' - User 14 bits per color component mode.
%
% Then we have support for the different modes of operation of the
% VPixx Technologies DPixx (DataPixx) box:
%
% 'M16' - Use 16 bit mono output mode, either with color index overlay
% (if the optional 2nd 'overlay' flag is set to 1, which is the default),
% or without color index overlay.
%
% 'C48' - User 16 bits per color component mode.
%
%
% 'BrightSide' - Enable drivers for BrightSide's HDR display. This only
% works if you have a BrightSide HDR display + the proper driver libraries
% installed on MS-Windows. On other operating systems it just uses a simple
% dummy emulation of the display with less than spectacular results.
%
% 'DualPipeHDR' - Use experimental output to dual-pipeline HDR display
% device.
%
%
% The third optional parameter 'colorclut' if provided, will use color
% lookup table based color correction / mapping, ie., mapping of luminance
% values to RGB values in the framebuffer - or intensity values if
% applicable. Any non-zero number will select a different CLUT. Look into
% the code for description. This just to demonstrate that you can use CLUT
% based color/intensity correction instead of simple power-law gamma
% correction if you want.
%
%
% The fourth optional parameter 'doGainCorrection' if provided and set to
% 1, will demonstrate use of display per-pixel gain correction, aka
% vignetting correction. It will modulate the brightness of each pixel with
% a gain factor, the gains increasing linearly from the left border to the
% right border of the display. See "help VignettingCorrectionDemo" for more
% details of this feature.
%
%
% Please note: Most of these modes only show expected results when the
% proper devices are attached and calibrated. All modes will work even on
% standard graphics without special devices, but you'll just see a false
% color image, as the standard GPU's can't interpret the special
% framebuffer encodings.
%
% The demo shows two superimposed sine wave gratings in the center of the
% screen. You can shift the 2nd grating up and down in subpixel steps by
% use of the cursor up-/down keys. You can change the contrast of the 2nd
% grating by use of the cursor left-/right keys. You can move the 2nd
% grating with the mouse, and rotate it clockwise or counterclockwise by
% pressing the mouse buttons. The keys 'i' and 'd' allow to change the
% "encoding gamma" factor used for the gamma correction algorithm. The ESC
% ape key ends the demo. Actually the demo performs some benchmark run for
% a few seconds after you've pressed ESC key, just to measure the speed of
% your graphics card in stimulus conversion.
%
% In 'VideoSwitcher' mode, it also draws some vertically moving greenish
% sync line just to show how to generate trigger signals on the
% VideoSwitcher device.
%
% 
% Needs hardware with support for imaging pipeline (GLSL shaders and
% floating point framebuffers). Should work well on ATI Radeon X1000 and
% later, Geforce 6000 and later and even better on DirectX10 hardware like
% Radeon HD series and NVidia Geforce 8 / 9 series and later.

% History:
% 16.04.2007 Written (MK).
% 13.08.2008 Cleaned up, commented, more help text etc... (MK).

KbName('UnifyKeyNames');
UpArrow = KbName('UpArrow');
DownArrow = KbName('DownArrow');
LeftArrow = KbName('LeftArrow');
RightArrow = KbName('RightArrow');
esc = KbName('ESCAPE');
space = KbName('space');
GammaIncrease = KbName('i');
GammaDecrease = KbName('d');

if nargin < 1 || isempty(outputdevice)
    outputdevice = 'None';
end

if nargin < 2 || isempty(overlay)
    overlay = 1;
end

if nargin < 3 || isempty(colorclut)
    colorclut = 0;
end

if nargin < 4 || isempty(doGainCorrection)
    doGainCorrection = 0;
end


try
	% This script calls Psychtoolbox commands available only in OpenGL-based 
	% versions of the Psychtoolbox. The Psychtoolbox command AssertPsychOpenGL will issue
	% an error message if someone tries to execute this script on a computer without
	% an OpenGL Psychtoolbox
	AssertOpenGL;
	
	% Get the list of screens and choose the one with the highest screen number.
	% Screen 0 is, by definition, the display with the menu bar. Often when 
	% two monitors are connected the one without the menu bar is used as 
	% the stimulus display.  Chosing the display with the highest dislay number is 
	% a best guess about where you want the stimulus displayed.  
	screenNumber=max(Screen('Screens'));
    
    % Open a double-buffered fullscreen window with a gray (intensity =
    % 0.5) background and support for 16- or 32 bpc floating point framebuffers.
    PsychImaging('PrepareConfiguration');

    lrect = [];
    rrect = [];
    
    % This will try to get 32 bpc float precision if the hardware supports
    % simultaneous use of 32 bpc float and alpha-blending. Otherwise it
    % will use a 16 bpc floating point framebuffer for drawing and
    % alpha-blending, but a 32 bpc buffer for gamma correction and final
    % display. The effective stimulus precision is reduced from 23 bits to
    % about 11 bits when a 16 bpc float buffer must be used instead of a 32
    % bpc float buffer:
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

    switch outputdevice
        case {'Mono++'}
            if overlay
                % Use Mono++ mode with color index overlay support:
                PsychImaging('AddTask', 'General', 'EnableBits++Mono++OutputWithOverlay');
            else
                % Use Mono++ mode without color overlay:
                PsychImaging('AddTask', 'General', 'EnableBits++Mono++Output');
            end

        case {'Color++'}
            % Use Color++ mode: We select averaging between even/odd
            % pixels, aka mode 2:
            PsychImaging('AddTask', 'General', 'EnableBits++Color++Output', 2);
            overlay = 0;

        case {'M16'}
            if overlay
                % Use M16 mode of Datapixx with color index overlay support:
                PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
            else
                % Use M16 mode of Datapixx without color overlay:
                PsychImaging('AddTask', 'General', 'EnableDataPixxM16Output');
            end

        case {'C48'}
            % Use C48 mode of Datapixx: We select averaging between even/odd
            % pixels, aka mode 2:
            PsychImaging('AddTask', 'General', 'EnableDataPixxC48Output', 2);
            overlay = 0;

        case {'Attenuator'}
            % Use the standard Pelli & Zhang style attenuator driver. This
            % uses a simple 3 row (for the three color channels Red, Green,
            % Blue) by n slots lookup table to map wanted intensity values
            % to RGB triplets for driving the attenuator. Any number of
            % slots up to 2^16 is supported, for a max precision of 16 bits
            % luminance. As we don't have a calibrated table in this demo,
            % we simply load a 2048 slot table (11 bit precision) with
            % random values:
            PsychImaging('AddTask', 'General', 'EnableGenericHighPrecisionLuminanceOutput', uint8(rand(3, 2048) * 255));
            overlay = 0;

        case {'VideoSwitcher'}
            % Select simple opmode of VideoSwitcher, where only the btrr
            % blue-to-red ratio from the global configuration file is used
            % for calibrated output:
            PsychImaging('AddTask', 'General', 'EnableVideoSwitcherSimpleLuminanceOutput', [], 1);
            
            % Switch the device to high precision luminance mode:
            PsychVideoSwitcher('SwitchMode', screenNumber, 1);
            overlay = 0;

        case {'VideoSwitcherCalibrated'}
            % Again the videoswitcher, but in lookup-table calibrated mode,
            % where additionally to the BTRR, a lookup table is loaded:
            PsychImaging('AddTask', 'General', 'EnableVideoSwitcherCalibratedLuminanceOutput', [], [], 1);
            
            % Switch the device to high precision luminance mode:
            PsychVideoSwitcher('SwitchMode', screenNumber, 1);
            overlay = 0;
            
        case {'PseudoGray'}
            % Enable bitstealing aka PseudoGray shader:
            PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
            overlay = 0;
            
        case {'Native10Bit'}
            % Enable ATI GPU's 10 bit framebuffer under certain conditions
            % (see help for this file):
            PsychImaging('AddTask', 'General', 'EnableNative10BitFramebuffer');
            overlay = 0;
        case {'BrightSide'}
            % Enable drivers for BrightSide's HDR display:
            PsychImaging('AddTask', 'General', 'EnableBrightSideHDROutput');
            overlay = 0;
            
        case {'None'}
            % No high precision output, just the plain 8 bit framebuffer,
            % but with software gamma correction:
            PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
            overlay = 0;

        case {'NoneNoGamma'}
            % No high precision output, just the plain 8 bit framebuffer,
            % even without gamma correction:
            PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
            overlay = 0;

        case {'DualPipeHDR'}
            % Enable experimental dual display, dual pipeline HDR output:
            
            % Handle single-screen vs. dual-screen output:
            if length(Screen('Screens')) == 1
                lrect = [0 0 600 600];
                rrect = [601 0 1201 600];
            end
            
            % Request actual output mode:
            PsychImaging('AddTask', 'General', 'EnableDualPipeHDROutput', min(Screen('Screens')), rrect);
            overlay = 0;

        otherwise
            error('Unknown "outputdevice" provided.');
    end

    if doGainCorrection
        % Request per pixel 2D gain correction for display:
        PsychImaging('AddTask', 'AllViews', 'DisplayColorCorrection', 'GainMatrix');
    end
    
    % Do not use gamma correction in calibrated video switcher mode or no
    % gamma mode -- wouldn't make sense in either of these:
    if ~strcmp(outputdevice, 'NoneNoGamma')
        if ~colorclut
            % Choose method of color correction: 'SimpleGamma' is simple gamma
            % correction of monochrome stims via power-law, ie., Lout = Lin ^ gamma.
            PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
            doTheGamma =1;
        else
            % Use CLUT based color correction:
            PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'LookupTable');
            doTheGamma =0;
        end
    else
        doTheGamma = 0;
    end
    
    %PsychImaging('AddTask', 'General', 'InterleavedLineStereo', 0);
    
    % Finally open a window according to the specs given with above
    % PsychImaging calls, clear it to a background color of 0.5 aka 50%
    % luminance:
    [w, wRect]=PsychImaging('OpenWindow',screenNumber, 0.5, lrect);
    
    % Use of overlay in Bits++ box Mono++ mode or DPixx box M16 mode wanted?
    if overlay
        % Get overlay window handle: Drawing into this window will affect
        % the overlay:
        wo = PsychImaging('GetOverlayWindow', w);
        Screen('Preference','TextAntiAliasing', 0);
    else
        wo = 0;
    end
    
    % Calibrated conversion driver for VideoSwitcher in use?
    if strcmp(outputdevice, 'VideoSwitcherCalibrated')
        % Tell the driver what luminance the background has. This allows
        % for some quite significant speedups in stimulus conversion:
        PsychVideoSwitcher('SetBackgroundLuminanceHint', w, 0.5);
    end

    % Assign a CLUT for color correction?
    if colorclut
        % Yes. 
        switch colorclut
            case 1,
                % Create a standard colormap() 3 channel RGB CLUT for a nice effect:
                clut = colormap;
            case 2,
                % A 1024 slots CLUT with inverted green channel:
                clut = ((0:1/1023:1)' * ones(1, 3));
                clut(:, 2) = 1 - clut(:, 2);
            case 3,
                % A simple inverted one channel (luminance) map with only two slots
                % for linearly interpolated values inbetween:
                clut = 1 - (0:1)';
            case 4,
                % Extreme amplification CLUT: Only meaningful on the
                % BrightSide HDR:
                clut = 3800 * (0:1)';
            otherwise,
                clut = colormap;
        end                

        % Assign the CLUT:
        PsychColorCorrection('SetLookupTable', w, clut);
    end
    
    if doTheGamma
        % We set initial encoding gamma to correct for a display with a
        % decoding gamma of 2.0 -- A good tradeoff, given most displays are
        % somewhere between 1.8 and 2.2:
        gamma = 1 / 2.0;
        PsychColorCorrection('SetEncodingGamma', w, gamma);
    else
        gamma = 1;
    end
        
    % From here on, all color values should be specified in the range 0.0
    % to 1.0 for displayable luminance values. Values outside that range
    % are allowed as intermediate results, but the final stimulus image
    % should be in range 0-1, otherwise result will be undefined.
    
    [width, height]=Screen('WindowSize', w);

    if doGainCorrection
        % Test support for display gain correction:
        gainmatrix = meshgrid(1:width, 1:height) / width * 1.5;
        PsychColorCorrection('SetGainMatrix', w, gainmatrix);
    end
    
    % Enable alpha blending. We switch it into additive mode which takes
    % source alpha into account:
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE);

	inc=0.5;
	
	% Compute one frame of a static grating: It has a total size of third
	% the screen size:
    s=min(width, height) / 6;
	[x,y]=meshgrid(-s:s-1, -s:s-1);
	angle=30*pi/180; % 30 deg orientation.
	f=0.01*2*pi; % cycles/pixel
    a=cos(angle)*f;
	b=sin(angle)*f;
                
    % Build grating texture:
    m=sin(a*x+b*y);
    tex=Screen('MakeTexture', w, m,[],[], 2);
    
    % Show the gray background:
    Screen('Flip', w);

    i=0;
    rotate = 0;
    yd =0;
    show2nd = 1;
    
    if overlay
        Screen('TextSize', wo, 18);
    else
        Screen('TextSize', w, 18);
    end
    
    % Center mouse on stimulus display, then make mouse cursor invisible:
    [cx, cy] = RectCenter(wRect);
    SetMouse(cx, cy, w);
    HideCursor;
    framecount = 0;
    
    tstart = GetSecs;
    
    % Animation loop:
    while 1
        Screen('DrawTexture', w, tex, [], [], [], [], 0.5);
        i=i+rotate;
        [x,y,buttons]=GetMouse(w);
        [x,y] = RemapMouse(w, 'AllViews', x, y);
        if any(buttons)
            if buttons(1)
                i=i+0.1;
            end
            if buttons(2)
                i=i-0.1;
            end
        end

        if show2nd
           dstRect=CenterRectOnPoint(Screen('Rect', tex), x, y+yd);
           Screen('DrawTexture', w, tex, [], dstRect, i, [], inc);
        end
        
        [d1 d2 keycode]=KbCheck;
        if d1
            if keycode(UpArrow)
                yd=yd-0.1;
            end

            if keycode(DownArrow)
                yd=yd+0.1;
            end

            if keycode(LeftArrow) && inc >= 0.001
                inc=inc-0.001;
            end

            if keycode(RightArrow) && inc <= 0.999
                inc=inc+0.001;
            end

            % Change of encoding gamma?
            if keycode(GammaIncrease) && doTheGamma
                gamma = min(gamma+0.001, 1.0);
                PsychColorCorrection('SetEncodingGamma', w, gamma);
            end
            
            if keycode(GammaDecrease) && doTheGamma
                gamma = max(gamma-0.001, 0.0);
                PsychColorCorrection('SetEncodingGamma', w, gamma);
            end

            if keycode(space);
                show2nd = 1-show2nd;
                while KbCheck; end;
                if show2nd
                    HideCursor;
                else
                    ShowCursor;
                end
            end

            if keycode(esc);
                break;
            end
        end

        txt0= 'At startup:\ngrating = sin(f*cos(angle)*x + f*sin(angle)*y);         % Compute luminance grating matrix in Matlab.\n';
        txt1= 'tex = Screen(''MakeTexture'', win, grating, [], [], 2); % Convert it into a 32bpc floating point texture.\n';
        txt2= 'Screen(''BlendFunction'', win, GL_SRC_ALPHA, GL_ONE);   % Enable additive alpha-blending.\n\nIn Display loop:\n\n';
        txt3= 'Screen(''DrawTexture'', win, tex, [], [], [], [], 0.5); % Draw static grating at center of screen.\n';
        txt4 = sprintf('Screen(''DrawTexture'', win, tex, [], [%i %i %i %i], %f, [], %f);\n', dstRect(1), dstRect(2), dstRect(3), dstRect(4), i, inc);
        txt5 = sprintf('\nEncoding Gamma is %f --> Correction for a %f gamma display.', gamma, 1/gamma);

        if overlay
            % Need to manually clear the overlay window:
            Screen('FillRect', wo, 0);
            % Need to use text color values in 0-255 range, instead of
            % normalized 0-1 range:
            DrawFormattedText(wo, [txt0 txt1 txt2 txt3 txt4 txt5], 0, 0, 255);
        else
            DrawFormattedText(w, [txt0 txt1 txt2 txt3 txt4 txt5], 0, 0, 1.0);
        end
        
        framecount = framecount + 1;
        
        % For the fun of it: Set a specific scanline to send a trigger
        % signal for the VideoSwitcher. This does nothing if the driver for
        % VideoSwitcher is not selected. We send out a trigger for 1 redraw
        % cycle every 30 redraw cycles. The triggerline is placed at
        % scanline 10 (for no special reason):
        if mod(framecount, 30) == 0
            PsychVideoSwitcher('SetTrigger', w, 10, 1);
        end
        
        % Show stimulus at next display retrace:
        Screen('Flip', w);
    end
        
    % Done.
    avgfps = framecount / (GetSecs - tstart);
    fprintf('Average redraw rate in demo was %f Hz.\n', avgfps);
    
    % Again, just to test conversion speed: A fast benchmark with sync of
    % buffer swaps to retrace disabled -- Go as fast as you can!
    nmaxbench = 300;
    tstart = Screen('Flip', w);
    for i=1:nmaxbench
        Screen('Flip', w, 0, 2, 2);
    end
    tend = Screen('Flip', w);
    fprintf('Average update rate in pipeline was %f Hz.\n', nmaxbench / (tend - tstart));
    
    Screen('Preference','TextAntiAliasing', 1);
    
    % We're done: Close all windows and textures:
    Screen('CloseAll');    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Screen('Preference','TextAntiAliasing', 1);
    ShowCursor;
    psychrethrow(psychlasterror);
end %try..catch..

if ~isempty(findstr(outputdevice, 'VideoSwitcher'))
    % If VideoSwitcher was active, switch it back to standard RGB desktop
    % display mode:
    PsychVideoSwitcher('SwitchMode', screenNumber, 0);
end

% Restore gfx gammatables if needed:
RestoreCluts;

return;
