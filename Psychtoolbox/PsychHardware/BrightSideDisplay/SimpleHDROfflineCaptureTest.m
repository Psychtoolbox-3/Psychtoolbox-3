function SimpleHDROfflineCaptureTest(imfilename, dummymode, sf)
% SimpleHDROfflineCaptureTest([imfilename][, dummymode][, sf])
%
% 'imfilename' - Filename of the HDR image to load. Will load our standard
% low-dynamic range 'konijntjes' if 'imfilename' is omitted.
%

% Make sure we run on OpenGL-Psychtoolbox. Abort otherwise.
PsychDefaultSetup(1);

% Run demo in dummy mode?
if nargin < 2
    dummymode = 0;
    if ~IsWin
        dummymode = 1;
    end
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
        
    imgRaw = [];
    
    % First render loop: Render each frame "live", then take a "screenshot"
    % of the raw HDR framebuffer image into a imgRaw() matrix for later
    % replay:
    
    % Initial Flip to sync us to retrace:
    vbl = Screen('Flip', win);
    tstart = vbl;
    
    for i=1:20
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
        rotAngle = rotAngle + 1;
        
        % Count our frames...
        framecounter = framecounter + 1;
        
        % Take raw HDR snapshot, return it as Matlab uint8 RGB matrix:
        imgRaw{framecounter} = BrightSideHDR('GetRawFramebufferSnapshot');

        %mcls = class(imgRaw{framecounter})
        %msiz = size(imgRaw{framecounter})
    end
      
    % We're done. Print the stats:
    framecounter
    duration = vbl - tstart
    averagefps = framecounter / duration

    Screen('Flip', win, vbl);
    
    rotAngle = 0;    
    framecounter = 0;
    imgDst = {};
    
    % Offline rendering: Create offscreen windows as buffers, render to
    % them and convert to HDR framebuffer formatted image, retrieve as
    % matrix for later replay:
    
    % Generate offscreen buffer pair for offscreen conversion of HDR
    % images:
    [sourceWin, destWin] = BrightSideHDR('CreateSnapshotBufferPair');
    
    
    % Now the same as offline render:
    for i=1:20
        % Clear backbuffer by overdrawing with a black full screen rect:
        Screen('FillRect', sourceWin, 0);

        % Draw our HDR texture at spec'd rotation angle:
        Screen('DrawTexture', sourceWin, texid, [], [], rotAngle);

        % Draw some 2D primitives:
        Screen('FillOval', sourceWin, [255 255 0], [500 500 600 600]);
        
        % And some text:
        Screen('TextSize', sourceWin, 30);
        Screen('TextStyle', sourceWin , 1);
        DrawFormattedText(sourceWin, 'Offline render demo.', 'center', 'center', [0 255 0]);
        
        % Increase rotation angle to make it a bit more interesting...
        rotAngle = rotAngle + 1;
        
        % Count our frames...
        framecounter = framecounter + 1;
        
        % Take raw HDR snapshot, return it as Matlab uint8 RGB matrix:
        imgDst{framecounter} = BrightSideHDR('ConvertImageToSnapshotBuffer', destWin, sourceWin);
        %[newfb, imgDst{framecounter}] = BrightSideHDR('GetSnapshotBuffer', sourceWin);

        %mcls = class(imgDst{framecounter})
        %msiz = size(imgDst{framecounter})
        
        % Playback immediately in real framebuffer:
        % BUG HERE: This should show the rendered image live, but instead
        % it shows a green or white screen. Unknown cause...
        % BrightSideHDR('BlitRawFramebufferSnapshot', destWin);
        % Screen('Flip', win);
        %BrightSideHDR('BlitRawFramebufferSnapshot', newfb);
        %Screen('Flip', win);
        %Screen('Close', newfb);
    end
    
    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');
    clear Screen;
    clear BrightSideHDR;
    
    % ================= SNIP ======================
    % Here you could end a session, storing the fetched imgRaw or imgDst
    % matrices to filesystem.
    %
    % Then later in a different script, load them back from filesystem,
    % then continue with code below for fast replay of such "prerecorded"
    % framebuffers:
    % ================= SNIP ======================
    
    
    
    % Try raw window + raw playback mode:
    win = BrightSideHDR('RawOpenWindow', screenid, 0);
    tic
    % Build raw HDR framebuffer textures from stored image matrices:
    for i=1:framecounter
        rawFB(i) = Screen('MakeTexture', win, imgRaw{i});
    end
    
    for i=1:framecounter
        rawFB2(i) = Screen('MakeTexture', win, imgDst{i});
    end
foo=    toc
    % Fast raw playback of first sequence of raw framebuffer textures:
    i = 0;
    while ~KbCheck
        i = mod(i+1, framecounter);
        
        BrightSideHDR('BlitRawFramebufferSnapshot', rawFB(i+1));
        Screen('Flip', win);        
    end
    
    KbReleaseWait;
    
    % Fast raw playback of 2nd sequence of raw framebuffer textures:
    i = 0;
    while ~KbCheck
        i = mod(i+1, framecounter);
        
        BrightSideHDR('BlitRawFramebufferSnapshot', rawFB2(i+1));
        Screen('Flip', win);        
    end
    
    % Close all textures, just to avoid superfluous warnings...
    Screen('Close');
    
    % Close everything else.
    Screen('CloseAll');
    
    % Restore gamma tables on HDR display:
    RestoreCluts;
    
    % Well done!
    fprintf('Bye bye!\n');
catch
    % Error handler: If something goes wrong between try and catch, we
    % close the window and abort.

    % Release all textures, close all windows, shutdown BrightSide library:
    Screen('CloseAll');

    % Restore gamma tables on HDR display:
    RestoreCluts

    % Rethrow the error, so higher level routines can handle it:
    psychrethrow(psychlasterror);
end;
