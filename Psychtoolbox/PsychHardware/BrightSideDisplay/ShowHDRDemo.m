function ShowHDRDemo(imfilename, dummymode, sf, halffloat)
% ShowHDRDemo([imfilename][, dummymode][, sf]) -- Load and show a high dynamic range image
% on the BrightSide Technologies High Dynamic Range display device.
%
% 'imfilename' - Filename of the HDR image to load. Will load our standard
% low-dynamic range 'konijntjes' if 'imfilename' is omitted.
%
% 'dummymode' - If set to 1 we only run in emulation mode without use of
% the HDR device or BrightSide core library.
%
% 'sf' - Scaling factor to apply.
%

% History:
% Written 2006 by Mario Kleiner - MPI for Biological Cybernetics, Tuebingen, Germany
% and Oguz Ahmet Akyuz - Department of Computer Science, University of Central Florida.

% Make sure we run on OpenGL-Psychtoolbox. Abort otherwise.
PsychDefaultSetup(1);

% Use remapping shaders? 0 = No, 1 = Yes. Shaders are no longer needed for
% HDR drawing, as we disable color clamping in the whole pipeline, even if
% in fixed-function mode.
useshader = 0;

% For now we shorten the sync tests during debugging:
Screen('Preference', 'SkipSyncTests', 1);

% Run demo in dummy mode?
if nargin < 2
    dummymode = 0;
end

% Name of an image file passed? If so, then load it. Load our default LDR
% image otherwise.
if nargin>=1 && ~isempty(imfilename)
    if ~isempty(strfind(imfilename, '.hdr'))
        % Load a real EXR high dynamic range file:
        inimg = read_rle_rgbe(imfilename);
        img(:,:,1) = flipud(inimg(:,:,1));
        img(:,:,2) = flipud(inimg(:,:,2));
        img(:,:,3) = flipud(inimg(:,:,3));
    else
        % Load a standard low dynamic range file:
        img = imread(imfilename);
    end
else
    % No. Just load our default low dynamic range image file:
    img = imread([PsychtoolboxRoot '/PsychDemos/konijntjes1024x768.jpg']);
end

% Scale by some factor:
if nargin < 3 || isempty(sf)
    if ~dummymode
        sf = 50;
    else
        sf = 1;
    end
end

% Enable half-float precision by default:
if nargin < 4 || isempty(halffloat)
    halffloat = 1;
end

% Scale image matrix by scalefactor 'sf':
img=double(img) * sf;

% img is now a height by width by 4 matrix with the (R,G,B,A) channels of our
% image. Now we use Psychtoolbox to make a HDR texture out of it and show
% it.
try
    % Find screen to display: We choose the one with the highest number,
    % assuming this is the HDR display:
    screenid=max(Screen('Screens'));

    % Initialize OpenGL mode of Psychtoolbox:
    InitializeMatlabOpenGL;
    
    % Open a standard fullscreen onscreen window, double-buffered with black
    % background color, instead of the default white one: win is the window
    % handle for this window:
    if dummymode
        win = BrightSideHDR('DummyOpenWindow', screenid, 0);
    else
        win = BrightSideHDR('OpenWindow', screenid, 0);
    end
    
    % Build a Psychtoolbox 16 bpc half-float or 32 bpc float texture from the image array:
    texid = Screen('MakeTexture', win, img, [], [], 2-halffloat);
    
    % Load our bias and rescale shader:
    glslnormalizer = LoadGLSLProgramFromFiles('ScaleAndBiasShader');
    prebias = glGetUniformLocation(glslnormalizer, 'prescaleoffset');
    postbias = glGetUniformLocation(glslnormalizer, 'postscaleoffset');
    scalefactor = glGetUniformLocation(glslnormalizer, 'scalefactor');

    % Activate it for setup:
    glUseProgram(glslnormalizer);

    % Set no bias to be applied:
    glUniform1f(prebias, 0.0);
    glUniform1f(postbias, 0.0);

    % Multiply all image values by a scaling factor: The HDR accepts
    % values between zero and infinity, but the useable range seems to be
    % zero (Dark) to something around 3000-4000. At higher values, it
    % saturates in a non-linear fashion. Fractional values, e.g,. 0.5 are
    % resolved at an unknown quantization level, so the range of
    % displayable intensity levels is more than 3000-4000 steps. We don't
    % know the real resolution without proper calibration, but according to
    % their Siggraph 2004 paper its supposed to be more than 14000 levels,
    % covering the full operating range in steps of single JND's. Who
    % knows... ... Ten is a good value for displaying LDR images...
    glUniform1f(scalefactor, 10);

    % Disable it. Will be enabled when needed:
    glUseProgram(0);

    % Load shader for normal primitive drawing (i.e. non-textures):
    glslcolor = LoadGLSLProgramFromFiles('HDRColorsShader');
    
    % Animation loop: Show a rotating HDR image, until user presses any key
    % to abort:
    rotAngle = 0;
    framecounter = 0;
    
    % Initial Flip to sync us to retrace:
    vbl = Screen('Flip', win);
    tstart = vbl;
    
    while ~KbCheck
        % Clear backbuffer by overdrawing with a black full screen rect:
        Screen('FillRect', win, 0);

        % Shall we use GLSL shaders?
        if useshader, glUseProgram(glslnormalizer); end;

        % Draw our texture into the backbuffer. We explicitely disable bilinear
        % filtering if halfffloat is 0, i.e., the texture is float32 format,
        % because current Geforce 7000 series hardware is not capable of
        % bilinear filtering of floating point 32 textures in hardware.
        % Bilinear filtering would work, but framerate would drop from 30
        % fps to 0.5 fps when the driver switches to the slow software
        % fallback path. If otoh halffloat is 1, i.e. the texture is
        % stored as 16 bpc half-floats, then we enable bilinear filtering,
        % as modern hardware is capable of filtering that.
        Screen('DrawTexture', win, texid, [], [], rotAngle, halffloat);

        % Draw some 2D primitives:
        if useshader
            glUseProgram(glslcolor);
        end
        
        Screen('FillOval', win, [255 * 255 * 10 255 0], [500 500 600 600]);
        
        % And some text:
        Screen('TextSize', win, 30);
        Screen('TextStyle', win , 1);
        if dummymode==0
            DrawFormattedText(win, 'If it works, it works.\nIf it doesn''t, it doesn''t.\n(Quoc Vuong, 2006)', 'center', 'center', [0 255*255*10 0]);
        end
        
        % Show updated HDR framebuffer at next vertical retrace:
        vbl=Screen('Flip', win, vbl);

        % Increase rotation angle to make it a bit more interesting...
        rotAngle = rotAngle + 0.1;
        
        % Count our frames...
        framecounter = framecounter + 1;
    end
  
    % We're done. Print the stats:
    framecounter
    duration = vbl - tstart
    averagefps = framecounter / duration
    
    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');
    
    % Well done!
    fprintf('Bye bye!\n');
catch
    % Error handler: If something goes wrong between try and catch, we
    % close the window and abort.

    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');

    % Rethrow the error, so higher level routines can handle it:
    psychrethrow(psychlasterror);
end;
