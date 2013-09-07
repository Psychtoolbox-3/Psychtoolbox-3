function SimpleHDRDemo(imfilename, dummymode, sf)
% SimpleHDRDemo([imfilename][, dummymode][, sf]) -- Load and show a high dynamic range image
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
% Btw. you can find lot's of nice free HDR images by Googling on the
% internet!
%
% Except for buying and installing the display device and control
% libraries, usage with Psychtoolbox is pretty straightforward. Modify your
% scripts in the following manner:
%
% 1. Use BrightSideHDR('OpenWindow', ...) instead of Screen('OpenWindow',
% ...) -- Will perform all additional display setup work for you.
%
% 2. Use HDRRead(imfilename) instead of imread(imfilename) to load HDR
% image files as Matlab matrices.
%
% 3. Set the 'floatprecision' flag of Screen('MakeTexture', ...) to 1 or 2
% to enforce creation of HDR textures from your image matrix.
%
% 4. Optionally you can use the Screen('ColorRange', ...) command to
% upscale or downscale color values for normal 2D drawing commands. This
% won't affect drawing of textures.
%
% 5. Optionally you can use fragment shaders to perform on-the-fly image
% processing when drawing an image, e.g., gamma correction, scaling, color
% conversion, tone-mapping. See the more advanced demos on how to do this.
%
% History:
% Written 2006 by Mario Kleiner - MPI for Biological Cybernetics, Tuebingen, Germany
% and Oguz Ahmet Akyuz - Department of Computer Science, University of Central Florida.

% Make sure we run on OpenGL-Psychtoolbox. Abort otherwise.
PsychDefaultSetup(1);

% Run demo in dummy mode?
if nargin < 2
    dummymode = 0;
end

% Name of an image file passed? If so, then load it. Load our default LDR
% image otherwise.
if nargin>=1 && ~isempty(imfilename)
    % Try to read a HDR image file by one of PTB's loaders. This will
    % return an empty img if format unsupported or unknown.
    img = HDRRead(imfilename, 1);
    
    if isempty(img)
        % Load a standard low dynamic range file:
        img = double(imread(imfilename));
    end
else
    % No filename provided. Just load our default low dynamic range image file:
    img = double(imread([PsychtoolboxRoot '/PsychDemos/konijntjes1024x768.jpg']));
end

% Scale intensities by some factor. Here we set the default value:
if nargin < 3 || isempty(sf)
    if ~dummymode
        sf = 50;
    else
        sf = 0.005;
    end
end

% Scale image matrix by scalefactor 'sf':
img=img * sf;

% img is now a height by width by c matrix with up to four (R,G,B,A) channels of our
% image. Now we use Psychtoolbox to make a HDR texture out of it and show
% it.
try
    % Find screen to display: We choose the one with the highest number,
    % assuming this is the HDR display:
    screenid=max(Screen('Screens'));
    
    % Open a standard fullscreen onscreen window, double-buffered with black
    % background color, instead of the default white one: win is the window
    % handle for this window. We use the BrightSideHDR() command instead of
    % Screen(). It is a convenience wrapper around Screen, doing all the
    % additional setup work for the HDR display:
    if dummymode
        % Dummy mode: Don't run on real HDR display:
        win = BrightSideHDR('DummyOpenWindow', screenid, 0);
    else
        % HDR mode: Setup everything for HDR rendering:
        win = BrightSideHDR('OpenWindow', screenid, 0);
    end
    
    % Build a Psychtoolbox 16 bpc half-float texture from the image array
    % by setting the (optional) 'floatprecision' flag to 1. If you need
    % even more precision you can provide the value 2 instead of 1,
    % creating full 32 bpc float textures. These will take up twice the
    % amount of memory and bandwidth though and they can't be anti-aliased
    % via bilinear filtering during drawing on current hardware - unless
    % you are happy with a framerate of 0.5 fps.
    texid = Screen('MakeTexture', win, img, [], [], 2);
    
    rotAngle = 0;
    framecounter = 0;
    
    % Initial Flip to sync us to retrace:
    vbl = Screen('Flip', win);
    tstart = vbl;
    
    % Animation loop: Show a rotating HDR image, until key press:
    while ~KbCheck
        % Clear backbuffer by overdrawing with a black full screen rect:
        Screen('FillRect', win, 0);

        % Draw our HDR texture at spec'd rotation angle:
        Screen('DrawTexture', win, texid, [], [], rotAngle);

        % Draw some 2D primitives:
        Screen('FillOval', win, [255 255 0], [500 500 600 600]);
        
        % And some text:
        Screen('TextSize', win, 30);
        Screen('TextStyle', win , 1);
        DrawFormattedText(win, 'If it works, it works.\nIf it doesn''t, it doesn''t.\n(Quoc Vuong, 2006)', 'center', 'center', [0 255 0]);
        
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
