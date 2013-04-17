function [win, winRect] = BrightSideHDR(cmd, arg, dummy, varargin)
% BrightSideHDR(cmd [, arg1][, arg2]) -- Psychtoolbox interface to
% BrightSide Technologies High Dynamic Range display device.
%
% This function is used to set up and interface with the High Dynamic Range
% display of BrightSide Technologies. It is a Matlab wrapper around lower
% level GLSL Psychtoolbox functions and the low-level BrightSideCore.dll
% MEX-File. You need to have Brightsides DLL's (and the display device)
% installed on your machine for this to work. Without the DLL's this will
% simply crash. See "help BSRuntimeLibs" for information about
% dependencies.
%
% cmd - The command that BrightSideHDR should execute. cmd can be any of
% the following:
%
% Open a window on the HDR display as with Screen('OpenWindow', ...),
% perform all HDR initialization:
%
% [win, winRect] = BrightSideHDR('OpenWindow', screenid, ...);
%
% This will execute Screen('OpenWindow') with all proper parameters,
% followed by BrightSide HDR init routines. It is a completely sufficient
% drop in replacement for Screen('OpenWindow'), accepting and returning
% exactly the same arguments that Screen() would do, adjusting all
% parameters to the constraints of the BrightSideHDR, if necessary.
%
% [win, winRect] = BrightSideHDR('DummyOpenWindow', screenid, ...);
% Same as OpenWindow, but just init a dummy mode that operates on standard
% displays without the Brightside-Libraries - useful for off-site testing.
%
%
% [win, winRect] = BrightSideHDR('RawOpenWindow', screenid, ...);
% Same as DummyOpenWindow, but just init a bare-minimum window on the HDR
% display suitable for playback of previously computed raw framebuffer
% images (see raw capture and replay support below). This will essentially
% open a fullscreen window on the HDR display, with proper window settings
% and identity gamma tables to make sure framebuffer content is passed
% unaltered 1-to-1 to the display. The imaging pipeline is disabled, this
% is a pure "fast replay" mode of content created in a previous session.
%
%
% In most applications you won't ever need to use any of the commands
% below, as PTB does all this stuff automatically behind the scenes...
% Ignore them unless you know what you're doing.
%
% BrightSideHDR('Initialize', win [, dummy]); -- Initialize the BrightSide libraries
% for HDR output into onscreen window 'win'. Onscreen window 'win' must have been
% opened before on a screen which corresponds to the attached High Dynamic Range
% display device. If you set 'dummy' to 1, then we are in emulation mode,
% i.e., we work without invocation of the mex file and without a real HDR
% display. This call will also attach all callback functions to the proper
% hooks for the given onscreen window 'win' inside the Screen() command.
%
% BrightSideHDR('Debuglevel', level); -- Set level of verbosity for
% debugging. The default is zero which means to be silent. A level of 1
% produces some debug output.
%
%
% SUPPORT FOR RAW HDR IMAGE CAPTURE AND REPLAY:
%
% [sourceWin, destWin] = BrightSideHDR('CreateSnapshotBufferPair');
% - Create a pair of offscreen windows useable for "offline" HDR->BrightSide
% conversion. Will create 'sourceWin' as a offscreen window with 32bpc
% float RGBA format, 'destWin' as a offscreenwindow with 8 bpc RGBA format.
% Both windows will have the size of the HDR displays framebuffer.
%
% You can draw into the 'sourceWin' just as you'd do to an onscreen window.
% Then you can call BrightSideHDR('ConvertImageToSnapshotBuffer', destWin,
% sourceWin); to convert the HDR image into a BrightSide formatted image in
% 'destWin'. You can use 'destWin' or its content in the
% BrightSideHDR('BlitRawFramebufferSnapshot') call for quick "replay" of
% HDR content.
%
%
% [rawImg] = BrightSideHDR('ConvertImageToSnapshotBuffer', destWin, sourceWin);
% - Convert HDR image in offscreen window 'sourceWin' into a HDR
% useable raw image in offscreen window 'destWin'. Optionally return a
% Matlab uint8 RGB image matrix with the raw data as 1st return argument
% rawImg.
%
% You could save 'rawImg' to disk and later load it and "replay" it on the
% HDR display via calls to BrightSideHDR('BlitRawFramebufferSnapshot',
% sourceWin);
%
%
% rawSnapshotMatrix = BrightSideHDR('GetRawFramebufferSnapshot');
% -Call immediately after Screen('Flip'). This will create a "screenshot" of the
% framebuffer for the currently displaying HDR image and return a Matlab
% matrix with proper raw framebuffer data, formatted for the HDR display.
% Can be used to Screen('MakeTexture') a texture that you can blit into the
% framebuffer via BrightSideHDR('BlitRawFramebufferSnapshot') without need
% for expensive conversion at runtime.
%
% 
% BrightSideHDR('BlitRawFramebufferSnapshot', sourceWin); Blit the given
% RGB8 or RGBA8 texture or offscreen window into the framebuffer of HDR
% display, so it gets displayed at the next invocation of Screen('Flip').
% The 'sourceWin' must be a texture or offscreen window whose image content
% was previously created by one of the raw snapshotting functions above,
% ie., it must be already encoded in the special data format for BrightSide
% displays.
%
%
% Callbacks: Usually called automatically by Screen(), no need to use them
% in your script directly! Some of them may disappear in the future,
% getting subsumed by Screen()'s internal imaging pipeline.
%
% BrightSideHDR('BrightSideExecuteBlit', win); -- This will convert
% the HDR image content in the HDR backbuffer into the special data
% format needed by the HDR display hardware. This is called by PTB
% automatically, don't call it yourself!
%
% BrightSideHDR('Shutdown'); -- Shut down the HDR device and core library,
% switch PTB back to standard LDR output. This is called by PTB
% automatically, don't call it yourself!

% History:
% 10/30/2006 Initial prototype implementation. MK & Oguz Ahmet Akyuz (Dept. of
% Computer Science, University of Central Florida)
% 11/09/2006 Small fixes and improvements to make it really work. (MK + Oguz)
% 12/12/2006 Integrate initial hook-plugin support for PTB's imaging
% pipeline, so user code does not need to call these functions anymore.
% 12/13/2006 Built-in 'OpenWindow' for simplified setup...
% 06/24/2007 Rewritten to work with current PTB, which has a full
% implementation of the floating point imaging pipeline. (MK)
% 07/05/2007 Bugfixes to rewrite. Improvements. Now its really ready for
% use with the new imaging pipeline. (MK)
% 07/01/2009 Add functions for capturing and replaying HDR framebuffers. (MK)

global GL;
persistent windowPtr;
persistent winwidth;
persistent winheight;
persistent hdrtexid;
persistent fboid;
persistent online;
persistent dummymode;
persistent debuglevel;
persistent oldTextRenderer;
persistent calledfromPsychImaging;
persistent reassign;

% Cold start: Setup our variables to safe defaults.
if isempty(online)
    online = 0;
    dummymode = 0;
    debuglevel = 0;
    hdrtexid = 0;
    reassign = 1;
    calledfromPsychImaging = 0;
end

if nargin < 1 || isempty(cmd)
    error('BrightSideHDR: Missing command specification.');
end

% Command dispatcher:
if strcmpi(cmd, 'Debuglevel')
    if nargin < 1 || isempty(arg)
        error('BrightSideHDR: "Debuglevel" called without specifiying a new level.');
    end
    
    debuglevel = arg;
    
    if ~dummymode && online
        % Set new debuglevel in core:
        BrightSideCore(-1, debuglevel);
    end
    
    return;
end

if strcmpi(cmd, 'CalledFromPsychImaging')
    % Mark us as being called from PsychImaging:
    if nargin < 1 || isempty(arg)
        error('BrightSideHDR: "CalledFromPsychImaging" called without specifiying a new level.');
    end
    calledfromPsychImaging = arg;

    return;
end

if strcmpi(cmd, 'OpenWindow') || strcmpi(cmd, 'DummyOpenWindow') || strcmpi(cmd, 'RawOpenWindow') || strcmpi(cmd, 'Initialize')

    % OpenGL mode of Psychtoolbox already initialized?
    if isempty(GL)
        % Nope. Do it now, but don't enable full 3D support, but only
        % support for calling MOGL commands, not z-/stencilbuffers and
        % context isolation:
        InitializeMatlabOpenGL([], [], 1);
    end

    % Open onscreen window as well or just perform initialization?
    if strcmpi(cmd, 'OpenWindow') || strcmpi(cmd, 'DummyOpenWindow') || strcmpi(cmd, 'RawOpenWindow')
        % Execute the Screen('OpenWindow') command with proper flags, followed
        % by our own Initialization. Return values of 'OpenWindow'.
        % BrightSideHDR('OpenWindow', ...) is a drop-in replacement for
        % Screen('OpenWindow', ...):


        % Assign screen index:
        if nargin < 2 || isempty(arg) || ~isa(arg, 'double')
            error('BrightSideHDR: "OpenWindow" called without valid screen handle.');
        end
        screenid = arg;
        
        % Assign optional clear color:
        if nargin < 3
            clearcolor = [];
        else
            clearcolor = dummy;
        end
        
        % windowRect is always full screen -- Anything else would make the
        % HDR display fail.
        winRect = [];
        
        % pixelSize is also fixed:
        pixelSize = 32;
        
        % Same for double-buffering:
        numbuffers = 2;
        
        % stereomode we take, unless raw window requested:
        if nargin >= 7 && ~strcmpi(cmd, 'RawOpenWindow')
            stereomode = varargin{4};
        else
            stereomode = [];
        end
        
        % multiSample gets overriden for now... Would probably interfere
        % with HDR display controller:
        multiSample = 0;
        
        % Imaging mode we take - and combine it with our own requirements:
        if nargin >= 9
            imagingmode = varargin{6};
        else
            imagingmode = 0;
        end
        
        % On a raw "playback only" window, we disable the imaging pipeline
        % to avoid interference:
        if ~strcmpi(cmd, 'RawOpenWindow')
            % Non-Raw: Need pipeline...
            
            % We need at least fast backing store and support for final output
            % conversion:
            imagingmode = mor(imagingmode, kPsychNeedFastBackingStore, kPsychNeedOutputConversion);

            % We require 32 bit floating point framebuffers if user doesn't
            % explicitely request the lower resolution, but higher speed 16 bit
            % floating point framebuffers:
            if bitand(imagingmode, kPsychNeed16BPCFloat) == 0
                imagingmode = mor(imagingmode, kPsychNeed32BPCFloat);
            end
        else
            % Raw window: Force pipeline off:
            imagingmode = 0;            
        end
        
        % Open the window, pass all parameters (partially modified or
        % overriden), return Screen's return values:
        if nargin > 9
            [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingmode, varargin{7:end});
        else
            [win, winRect] = Screen('OpenWindow', screenid, clearcolor, winRect, pixelSize, numbuffers, stereomode, multiSample, imagingmode);
        end
        
        % Ok, if we reach this point then we've got a proper onscreen
        % window on the HDR. Let's reassign our arguments and continue with
        % the init sequence:
        arg = win;
        
        % Operate in dummy mode or real mode?
        if strcmpi(cmd, 'DummyOpenWindow') || strcmpi(cmd, 'RawOpenWindow')
            dummy = 1;
        else
            dummy = 0;
        end
        
        % Ok, should be safe to continue with standard Brightside init...
    end

    % BrightSide init sequence:

    % Child protection:
    if online
        error('BrightSideHDR: "Initialize" called on a HDR display already initialized.');
    end

    if nargin < 2 || isempty(arg) || ~isa(arg, 'double')
        error('BrightSideHDR: "Initialize" called without valid HDR onscreen window handle.');
    end

    if strcmpi(cmd, 'Initialize') && (nargin < 3 || isempty(dummy))
        dummymode = 0;
    else
        dummymode = dummy;
    end
    
    % Another bit of child protection:
    AssertGLSL;

    % Step 1: Make sure we got a window handle for a valid onscreen window:
    if isempty(find(Screen('Windows')==arg, 1)) || Screen('WindowKind', arg)~=1
        error('BrightSideHDR: "Initialize" called with something else than a valid onscreen window handle.');
    end
    
    % Valid handle, assign it and query its dimensions, aka the dimensions
    % of our HDR framebuffer:
    windowPtr = arg;
    [winwidth, winheight] = Screen('WindowSize', windowPtr);
        
    % Only setup the BrightSide drivers and libraries in non-dummy mode:
    if ~dummymode
        % Initiate loading, linking and initialization of the core:

        % On Windows, we need to preload BrightSideCore into Matlab while the working
        % directory is set to Psychtoolbox/PsychHardware/BrightSideDisplay/BSRuntimeLibs/outputlib/lib/
        % , so the Windows dynamic linker can find our own local copies of the BrightSide DLL's.
        % and link against it.
        if IsWin
            % Windows system: Change working dir to location to BrightSide
            % DLL's:
            olddir = pwd;
            
            try
                % This fails if libraries are not installed, so we
                % try-catch it...
                cd([PsychtoolboxRoot 'PsychHardware/BrightSideDisplay/BSRuntimeLibs/outputlib/lib']);
            catch
                fprintf('Could not cd() into the directory %s !\n', [PsychtoolboxRoot 'PsychHardware/BrightSideDisplay/BSRuntimeLibs/outputlib/lib']);
                fprintf('This could be a problem with access permissions, but more likely, this directory simply does not exist.\n');
                fprintf('After a new installation of Psychtoolbox, you must manually copy the BrightSide runtime libraries into that\n');
                fprintf('folder, ie., replace the empty directory .../BSRuntimeLibs/ by the one contained in the zip file with the\n');
                fprintf('proprietary Brightside runtime libraries.\n');
                fprintf('As of June 2008, in MPI Tuebingen, this zip file can be found in \\\\uni.kyb.local\\home\\kleinerm\\projects\\OpenGLPsychtoolbox\\BrightSideHDRRuntime.zip\n');
                fprintf('Check this, then retry...\n');
                
                Screen('CloseAll');
                error('BrightSideHDR: The directory with the Brightside runtime drivers is inaccessible or non-existent.');
            end
            
            % Preload (and thereby link) BrightSideCore into Matlab. The
            % special command code -1 forces loading and sets the initial
            % debuglevel:
            BrightSideCore(-1, debuglevel);

            % Now that it is (hopefully) properly loaded, we can revert the working
            % directory to its previous setting:
            cd(olddir);
        end;        

        % Startup & Initialize the BrightSideHDR driver library. We pass the
        % path to all configuration files and the name of the master
        % configuration file to the startup routine.
        BrightSideCore(0, [fileparts(which('BrightSideCore')) '/BSRuntimeLibs/Resources'], 'DR-37P-beta.xml');

        % Add the BrightSideHDR blitter callback to the final output
        % formatter blitchain: This will trigger data conversion at each
        % invocation of Screen('Flip'):
        Screen('HookFunction', windowPtr, 'AppendMFunction', 'FinalOutputFormattingBlit', 'Execute BrightSide blit operation', 'BrightSideHDR(''BrightSideExecuteBlit'')');
        Screen('HookFunction', windowPtr, 'Enable', 'FinalOutputFormattingBlit');
        
        % Load the gfx-hardware gamma-tables / CLUT's with identity, so
        % they can't mess with BrightSide encoded image data:
        LoadIdentityClut(windowPtr);
    end

    % Add a shutdown callback - Always called (even in dummy mode) to
    % teardown the BrightSide-Lib (if active) and to reset our internal
    % data structures:
    Screen('HookFunction', windowPtr, 'AppendMFunction', 'CloseOnscreenWindowPreGLShutdown', 'Shutdown BrightSide core before window close.', 'BrightSideHDR(''Shutdown'')');
    Screen('HookFunction', windowPtr, 'Enable', 'CloseOnscreenWindowPreGLShutdown');

    if ~strcmpi(cmd, 'RawOpenWindow')
        % Disable color clamping in the GL pipeline. It's not useful for our
        % purpose, unless we are in dummymode. Also set color scaling to 1.0,
        % i.e. do not scale color values at all - Doesn't make much sense with
        % the BrightSide HDR display.
        Screen('ColorRange', windowPtr, 1, 0);

        % Eat up all OpenGL errors caused by this:
        while glGetError; end;
    else
        % A raw window for pure playback of previously "recorded" buffers
        % still needs a linearized video output to work properly:

        % Load the gfx-hardware gamma-tables / CLUT's with identity, so
        % they can't mess with BrightSide encoded image data:
        LoadIdentityClut(windowPtr);
    end
    
    % Use old display list text renderer - The new one is not HDR ready: It
    % has a restricted color value range 0.0-1.0 for text colors, which
    % will create way too dark text, and its internal use of alpha-blending
    % triggers software-fallback in the Geforce 7000 gfx... Thereby, use
    % the old and trouble-free one:
    oldTextRenderer = Screen('Preference', 'TextRenderer', 0);
    
    % We are online:
    online = 1;

    % We must skip the following code if called from PsychImaging!
    % Reason: PsychImaging will execute this call sequence internally, but
    % after setup of its PsychColorCorrection - ie., in a different order.
    % The first call to Screen('Flip') will determine which internal FBO is
    % the correct source FBO for BrightSide libs, and PsychImaging may
    % still change/add/swizzle the order of internal FBO's and imaging
    % operations, so it would be too early for us to bind the final source
    % FBO, iow. to do the first 'Flip'. If we 'flipped' here, the wrong FBO
    % might be permanently bound, thereby preventing proper image
    % processing sequence in PsychImaging code:
    if ~calledfromPsychImaging
        % Set the background clear color via old fullscreen 'FillRect' trick,
        % followed by a flip:
        Screen('FillRect', windowPtr, clearcolor);
        Screen('Flip', win);
    end
    
    return;
end

if strcmpi(cmd, 'Shutdown')
    % Shutdown command:

    % Child protection:
    if ~online
        warning('BrightSideHDR: "Shutdown" command called, although display is already offline.'); %#ok<WNTAG>
        return;
    end
    
    % Shutdown BrightSide core and display:
    if ~dummymode
        BrightSideCore(1);
    end

    % Restore selection of text renderer:
    % TODO FIXME: Technically it is not allowed to call into Screen() from
    % within this function as it is called by Screen()'s shutdown path and
    % Screen is not reentrant! For some reason it still works. Are we lucky
    % or what?!? Leave it for now...
    Screen('Preference', 'TextRenderer', oldTextRenderer);

    % Reset to preinit state:
    hdrtexid = 0;
    fboid = 0;
    windowPtr = 0;
    dummymode = 0;
    reassign = 1;

    % Shutdown complete, we're offline:
    online = [];    
    return;
end

if strcmpi(cmd, 'BrightSideExecuteBlit')
    % Initiate float FBO -> Backbuffer data conversion to convert a HDR
    % float image into data format required by BrightSide HDR:
    if ~online
        error('BrightSideHDR: "BrightSideExecuteBlit" command called, although display is offline.');
    end

    % If we are still in the PsychImaging setup phase...
    if calledfromPsychImaging
        % Then turn this call into a no-op!
        return;
    end
    
    % Disable shaders, if any active:
    glUseProgram(0);

    % 'hdrtexid' of source color texture already assigned?
    if hdrtexid == 0
        % No. But the texture currently bound to unit 0 is the proper one,
        % so query its id. We know that the texture is bound to the
        % rectangle texture target and that unit 0 is active, due to the
        % semantics of PTB's imaging pipeline. This subroutine is called
        % from the 'OutputFormattingBlit'-Hookchain of Screen('Flip'). Part
        % of the hookchain execution sequence is setup of proper source
        % textures and target FBOs, so we know the currently bound texture
        % is the floating point source texture and the currently bound
        % framebuffer is the target framebuffer (FBO or system
        % framebuffer):
        hdrtexid = double(glGetIntegerv(GL.TEXTURE_BINDING_RECTANGLE_EXT));
        fboid = double(glGetIntegerv(GL.FRAMEBUFFER_BINDING_EXT));

        if hdrtexid <= 0
            error('Fatal error: Failed to query color buffer texture id in BrightSideHDR("BrightSideExecuteBlit")!');
        end
        
        if fboid < 0
            error('Fatal error: Failed to query Framebuffer object id in BrightSideHDR("BrightSideExecuteBlit")!');
        end

        % Need to do first-time init of HDR libary with this 'hdrtexid':
        needhdrinit = 1;
    else
        % No need to first-time init HDR library:
        needhdrinit = 0;
    end
    
    % Execute actual conversion op:
    DoExecuteHDRBlit(dummymode, winwidth, winheight, needhdrinit, hdrtexid, fboid);
    
    % Ready:
    return;
end

if strcmpi(cmd, 'CreateSnapshotBufferPair')
    if ~online
        error(sprintf('BrightSideHDR: "%s" command called, although display is offline.', cmd)); %#ok<SPERR>
    end

    % Create source offscreen window: Has 128 bpp or 32bpc float RGBA
    % representation:
    srcBuffer = Screen('OpenOffscreenWindow', windowPtr, 0, [], 128, 32);

    % Yes. Create a target RGBA8 offscreen window of sufficient size:
    dstBuffer = Screen('OpenOffscreenWindow', windowPtr, 0, [], 32, 32);
    
    win = srcBuffer;
    winRect = dstBuffer;
    
    return;
end

% Convert a HDR floating point image (inside a texture or offscreen window)
% into an RGBA8 raw image buffer, suitable for directly driving the HDR display.
if strcmpi(cmd, 'ConvertImageToSnapshotBuffer') || strcmpi(cmd, 'GetSnapshotBuffer')
    if ~online
        error(sprintf('BrightSideHDR: "%s" command called, although display is offline.', cmd)); %#ok<SPERR>
    end

    % Shall we create a suitable target RGBA8 buffer ourselves?
    if strcmpi(cmd, 'GetSnapshotBuffer')
        % Yes. Create a RGBA8 offscreen window of sufficient size:
        dstBuffer = Screen('OpenOffscreenWindow', windowPtr, 0, [], 32, 32);
        
        if exist('arg', 'var') && ~isempty(arg)
            srcBuffer = arg;            
        else
            % No srcBuffer provided:
            srcBuffer = [];
        end

        scmode = 1;
        
        % First return argument will be destination offscreen window handle:
        win = dstBuffer;        
    else
        % No. Sanity check and then use the provided handle:
        if ~exist('arg', 'var') || isempty(arg)
            error(sprintf('BrightSideHDR:%s: Required destinationWindow handle missing!', cmd)); %#ok<SPERR>
        end
        
        % Assign dstBuffer handle:
        dstBuffer = arg;
        
        % Proper type of destination?
        if Screen('WindowKind', dstBuffer) ~= -1
            error(sprintf('BrightSideHDR:%s: Invalid destinationWindow handle passed. Not an offscreen window!', cmd)); %#ok<SPERR>
        end

        % Matching dimensions with display framebuffer?
        [w, h] = Screen('WindowSize', dstBuffer);
        if w~=winwidth || h~=winheight
            error(sprintf('BrightSideHDR:%s: Size of destinationWindow %i x %i pixels doesn''t match size of display %i x %i pixels as required!', cmd, w, h, winwidth, winheight)); %#ok<SPERR>
        end
        
        % Check format of destination window: Must be a 8bpc, RGBA8, 32
        % bpp:
        bpp = Screen('PixelSize', dstBuffer);
        if  bpp ~= 32
            error(sprintf('BrightSideHDR:%s: Pixelsize %i bpp of destinationWindow handle is not the required 32 bpp for a RGBA8 target!', cmd, bpp)); %#ok<SPERR>
        end
        
        % Ok, checks passed. Check optional source argument:
        if exist('dummy', 'var') && ~isempty(dummy)
            % Assign it:
            srcBuffer = dummy;  
        else
            % No srcBuffer provided:
            srcBuffer = [];
        end
        scmode = 0; 
    end
    
    % Extract srcTex texture handle for the RGBA floating point source
    % texture for BrightSide conversion:
    
    % srcBuffer provided?
    if isempty(srcBuffer)
        % No source window provided: We default the srcTex texture id
        % to the texture of our regular PTB imaging pipeline
        % framebuffer:
        srcTex = hdrtexid;
        
        % Make sure pending rendering operations are finished:
        Screen('DrawingFinished', windowPtr, [], 1);
    else
        % Proper type of source?
        if Screen('WindowKind', srcBuffer) ~= -1
            error(sprintf('BrightSideHDR:%s: Invalid sourceWindow handle passed. Not an offscreen window or texture!', cmd)); %#ok<SPERR>
        end

        % Matching dimensions with display framebuffer?
        [w, h] = Screen('WindowSize', srcBuffer);
        if w~=winwidth || h~=winheight
            error(sprintf('BrightSideHDR:%s: Size of sourceWindow %i x %i pixels doesn''t match size of display %i x %i pixels as required!', cmd, w, h, winwidth, winheight)); %#ok<SPERR>
        end

        % srcBuffer is eligible as source window for conversion.
        % Extract its OpenGL texture id:
        srcTex = Screen('GetOpenGLTexture', windowPtr, srcBuffer);
        
        % Make sure pending rendering operations are finished:
        Screen('DrawingFinished', srcBuffer, [], 1); 
    end

    % Need to extract FBO id for dstBuffer of BrightSide conversion:
    
    % The call to 'GetWindowInfo' will as a side effect bind the FBO for dstBuffer:
    Screen('GetWindowInfo', dstBuffer);

    % Query it:
    dstFBOid = double(glGetIntegerv(GL.FRAMEBUFFER_BINDING_EXT));
    
    % Disable shaders, if any active:
    glUseProgram(0);
        
    if ~dummymode
        % Ok, we have a valid source texture id and destination FBO id. Call
        % BrightSideCore to setup the library for conversion between these buffers:
        DoExecuteHDRBlit(0, winwidth, winheight, reassign, srcTex, dstFBOid);
    else
        % In dummymode, DoExecuteHDRBlit() is a no-op, so we need to
        % replace it by a plain 1->1 blit:
        Screen('CopyWindow', srcBuffer, dstBuffer);
    end
    % Switch library back to standard source and destination assignment:
    if ~dummymode
        % MK BUG: reassign should be 1 and the commented out command should be
        % used: We should switch back and forth in the assignment of source
        % textures and target fbos for the BrightSide library for optimal
        % behaviour. However the library doesn't like multiple successive
        % switches at all -- The whole beast fails if we try to do it.
        % Therefore we perform this assignment only 1 time --> Doesn't
        % allow us to switch between offline and online use!
        % It is unclear if we encounter a bug in the BrightSide library,
        % some Psychtoolbox bug, some coding bug here or a NVidia driver
        % bug. We do know however that the BrightSide libs are often
        % leaking OpenGl state and are pretty bad and faulty in managing state, so it
        % might be a BrightSide bug...
        reassign = 0;
        %        BrightSideCore(3, hdrtexid, fboid);
    end
    
    % Ok, now we should have our final BrightSide ready RGB8 framebuffer
    % image in dstBuffer. Should we also return a snapshot as Matlab RGB8
    % matrix in uint8 format, suitable for later texture creation?
    if (scmode && nargout > 1) || (~scmode && nargout > 0)
        % Yes. Get "screenshot"
        outImg = Screen('GetImage', dstBuffer, [], [], 0, 3);
        
        % Assign to proper return argument:
        if scmode
            winRect = outImg;
        else
            win = outImg;
        end
    end

    % Done.
    return;
end

if strcmpi(cmd, 'GetRawFramebufferSnapshot')
    if ~online
        error(sprintf('BrightSideHDR: "%s" command called, although display is offline.', cmd)); %#ok<SPERR>
    end
    
    % Retrieve RGB8 uint8 snapshot image matrix of the system frontbuffer
    % suitable for creating a raw RGB8 texture for driving the BrightSide
    % display without imaging pipeline, kind'a raw playback:
    win = Screen('GetImage', windowPtr, [], 'frontBuffer', 0, 3);
    
    % Done.
    return;
end

if strcmpi(cmd, 'BlitRawFramebufferSnapshot')
    if ~online
        error(sprintf('BrightSideHDR: "%s" command called, although display is offline.', cmd)); %#ok<SPERR>
    end

    if ~exist('arg', 'var') || isempty(arg)
        error(sprintf('BrightSideHDR:%s: Mandatory srcWindow handle argument missing!', cmd)); %#ok<SPERR>
    end
    
    srcBuffer = arg;
    
    % Proper type of source?
    if Screen('WindowKind', srcBuffer) ~= -1
        error(sprintf('BrightSideHDR:%s: Invalid sourceWindow handle passed. Not an offscreen window or texture!', cmd)); %#ok<SPERR>
    end

    % Matching dimensions with display framebuffer?
    [w, h] = Screen('WindowSize', srcBuffer);
    if w~=winwidth || h~=winheight
        error(sprintf('BrightSideHDR:%s: Size of sourceWindow %i x %i pixels doesn''t match size of display %i x %i pixels as required!', cmd, w, h, winwidth, winheight)); %#ok<SPERR>
    end
    
    % Make sure that at least the final output formatter is off:
    Screen('HookFunction', windowPtr, 'Disable', 'FinalOutputFormattingBlit');
    
    % Perform blit: We use CopyWindow to get an identity blit, possibly
    % accelerated via FBO blit extensions on modern hardware:
    %Screen('CopyWindow', srcBuffer, windowPtr);
        Screen('DrawTexture', windowPtr, srcBuffer, [], [], [], 0);
    
    % We enforce finalization of preflip operations:
    Screen('DrawingFinished', windowPtr, [], 1);
        
    % Reenable output formatter if it was active before:
    if Screen('HookFunction', windowPtr, 'Query', 'FinalOutputFormattingBlit', 0) ~= -1
        Screen('HookFunction', windowPtr, 'Enable', 'FinalOutputFormattingBlit');
    end

    % Done.
    return;
end

% Unknown command.
error('BrightSideHDR: Unknown subcommand specified! Type ''help BrightSideHDR'' for usage info.');
end

% Helper routine: Performs actual library calls to perform HDR -> RGBA8 conversion:
function DoExecuteHDRBlit(dummymode, winwidth, winheight, needhdrinit, mhdrtexid, mfboid)
global GL;

    if ~dummymode
        % Setup special projection matrices for the BrightSide core lib:
        glMatrixMode(GL.PROJECTION);
        glPushMatrix;
        glLoadIdentity;

        % The following gluOrhto2D command causes an OpenGL error, because
        % it defines invalid settings, but for some reason we need exactly this
        % transformation for the BrightSideHDR to work. Therefore we disable
        % automatic error-checking in moglcore, then execute the command, then
        % eat up the error-code and then reenable error checking:

        % Disable automatic error checking:
        moglcore('DEBUGLEVEL', 0);

        % Perform the strange setup command:
        gluOrtho2D(0, 0, winwidth, winheight);

        % Eat up all OpenGL errors caused by this:
        while glGetError; end;

        % Reenable automatic error checking:
        moglcore('DEBUGLEVEL', 1);

        glMatrixMode (GL.MODELVIEW);
        glPushMatrix;
        glLoadIdentity;

        % Call the BrightSide core library to do the conversion:
        % It is crucial to backup & restore at least the GL_ENABLE_BIT
        % state around the call to BrightSides core libraries. They seem to
        % change enable-state and don't restore it properly. For that
        % reason, we couldn't use unclamped colors when the fixed-function
        % pipeline was active, i.e. no shaders bound.        
        glPushAttrib(GL.ALL_ATTRIB_BITS);

        % Would be enough to push this, but better safe than sorry. You never
        % know what your software vendor will screw in its next release ;)
        % glPushAttrib(GL.ENABLE_BIT);

        if needhdrinit == 1
            % Setup of source texture and target FBO for HDR -> Framebuffer
            % conversion. We do it here because we need the 'hdrtexid' of the
            % source color texture, which is only reliably available at first
            % invocation of this subroutine. Same for 'fboid' (although
            % this will be zero in 99.99% of all PTB applications ;-) )
            BrightSideCore(3, mhdrtexid, mfboid);
        end
        
        % Perform actual conversion-blit:
        BrightSideCore(2);

        % Restore saved state:
        glPopAttrib;

        % Restore normal matrices:
        glMatrixMode(GL.PROJECTION);
        glPopMatrix;
        glMatrixMode (GL.MODELVIEW);
        glPopMatrix;
    end
    
    return;
end
