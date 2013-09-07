function varargout = PsychVideoSwitcher(cmd, varargin)
% PsychVideoSwitcher(command [,arg1, arg2, .....]);
%
% Psychtoolbox support for the Xiangrui Li et al. "VideoSwitcher" video
% attenuator device for high precision luminance output with up to 16 bits
% luminance resolution.
%
% This routine incorporates code contributed by Xiangru Li for switching
% between monochrome and color display mode and for performing the
% reference Matlab routines for image formatting.
%
% Options: 'command' is a command string, specifying a subcommand. The
% following subcommands are supported with the following options:
%
% Public functions:
%
%
% PsychVideoSwitcher('SwitchMode', screenIdx, enableLuminanceMode [, VideoSwitcherIsABox])
% - Switch programmatically between high precision luminance mode and
% standard RGB true color display mode. 'screenIdx' is the screen index of the
% for the display screen to switch.
%
% 'enableLuminanceMode' is meaningful only for card VideoSwitcher, which
% must be set to 0 to switch to RGB mode, and to 1 to switch to high
% precision luminance mode. For box version, the SwitchMode subfunction
% ignores the third input, and only toggles between two display modes,
% equivalent to pushing the switch button on the box.
%
% 'VideoSwitcherIsABox' is an optional argument: If set to 1, then perform
% switching procedure for an external (box) device. If set to zero, then
% perform procedure for an internal (PCI card) device. If argument is
% omitted, the proper default type is read from a configuration file stored
% in the path mypath = PsychtoolboxConfigDir('VideoSwitcher');
%
% If you create a file named 'VideoSwitcherIsABox' in that directory
% mypath, the switcher is assumed to be an external box. If you create a
% file named 'VideoSwitcherIsACard', the code assumes a PCI card based
% device. If no such files exist, the switcher is assumed to be a box.
%
% It is important for the code to know if the switcher is a box or a card,
% as the switching strategy is different.
%
%
% PsychVideoSwitcher('SetTrigger', win, triggerLine [, count=infinite]);
% - Set trigger line and options for VideoSwitcher connected to the display
% of onscreen window 'win'.
%
% 'triggerLine' defines the vertical (y) position of a trigger line to be
% drawn to the green channel of the final image, in order to trigger the
% trigger-circuit of the VideoSwitcher. If you set it to a negative value
% or to empty [], triggering will be disabled (This is the default).
%
% 'count' Optional: Number of redraw cycles (invocations of
% Screen('Flip')), that the triggerline should be drawn. If left out, the
% line will be drawn on each redraw until disabled. If set to some value,
% it will be drawn for 'count' number of redraws, then disabled.
%
% The trigger mechanism only works if the VideoSwitcher is set to high
% precision luminance display. It will emit a TTL pulse of 100 microseconds
% duration when the triggerline is drawn by your display. Only one TTL
% trigger pulse per video refresh is possible.
%
%
% PsychVideoSwitcher('SetBackgroundLuminanceHint', win, luminance);
% - Tell the driver the 'luminance' value of the background pixels of
% onscreen windows 'win'. The driver will use this hint to optimize
% conversion of background pixels. This allows for a quite significant
% speedup if your stimulus only covers a fraction of the display area and
% the remaining area is just a uniform luminance background. This is only
% functional if you use the driver that uses the calibration LUT, not the
% simple driver.
%
%
% RGBImage = PsychVideoSwitcher('MapLuminanceToRGB', lum, ratio [, trigger]);
% - Perform conversion of a luminance image into a RGBImage. This is a pure
% Matlab based implementation for graphics hardware that is not capable of
% supporting the imaging pipeline.
%
%  Inputs:
%     lum: luminance (MxN matrix with values from 0 to 1)
%     ratio: blue to red ratio of the video switcher
%     trigger: when non-zero, a trigger will be sent in current frame
%         trigger=1 or 'top', the first line of image
%         trigger=2 or 'auto',  the first line with non-zero image
%         trigger=3 or 'middle', the middle line of image
%     If you omit ratio, you should give it in this code
%  Output: RGB image (MxNx3 matrix with values from 0 to 255)
%
%
% RGBImage = PsychVideoSwitcher('MapLuminanceToRGBCalibrated', lum, ratio, lut [, trigger])
% - Perform conversion of a luminance image into a RGBImage. This is a pure
% Matlab based implementation for graphics hardware that is not capable of
% supporting the imaging pipeline.
%
%  Inputs:
%     lum: luminance (MxN matrix with values from 0 to 1)
%     ratio: blue to red ratio of the video switcher
%     lut: The 257 slots calibrated luminance table as described below.
%     trigger: when non-zero, a trigger will be sent in current frame
%         trigger=1 or 'top', the first line of image
%         trigger=2 or 'auto',  the first line with non-zero image
%         trigger=3 or 'middle', the middle line of image
%     If you omit ratio, you should give it in this code
%  Output: RGB image (MxNx3 matrix with values from 0 to 255)
% 
% Note: this runs slow when lum matrix is large. 200x200 can take 1 second.
%
% This uses the calibrated luminance table saved in calibratedlum.mat.
% For detail, check the paper http://lobes.usc.edu/Journals/JNM03.pdf
% 
% The steps to get the table:
%
%  1. Switch the video switcher to grayscale mode;
%
%  2. Set up the equipement to accurately measure screen luminance;
%     you can use a photometer or data acquisition system;
%
%  3. Measure 257 luminance levels at RGB of [0 0 b] and [btrr 0 255],
%     where b is 0:255, and btrr the blue to red ratio of switcher; 
%
%  4. Store 257 luminance in a variable callum and normalize them:
%     callum=callum/callum(257);
%
%  5. Save it to a configuration file. See the help for subfunction
%     'GetDefaultConfig' below for how and where to store the calibration
%     table.
%
%
% [btrr, lut] = PsychVideoSwitcher('GetDefaultConfig', win);
% - Get default 'btrr' parameter and 'lut' lookup table from configuration
% files and return them. This function can be used by you to get switcher
% parameters for use with the Matlab conversion functions above. It will be
% automatically used by the imaging pipeline (by PsychImaging() command) if
% you select VideoSwitcher as an output device but don't provide explicit
% settings for 'btrr' and/or 'lut'.
%
% Configuration files for the VideoSwitcher should be stored in the folder
% whose path you get if you type mypath = PsychtoolboxConfigDir('VideoSwitcher');
%
% You can store a configuration file specific to a display screen
% 'screenid' (with the numbering as in Screen('Screens')). The file should
% have the name SettingsforScreen_X.mat with X being the screen number,
% e.g., SettingsforScreen_0.mat . This allows you to store per-display
% device settings. Alternatively you can store settings in a "global"
% config file named GlobalSettings.mat if they are not specific to the
% display. The stored mat file should contain up to two variables:
%
% The variable 'btrr' should store the measured/calibrated BTRR
% Blue-To-Red-Ratio for yor switcher and display setup. This variable is
% mandatory.
%
% The optional variable 'lut' would be a 257 elements double vector of
% luminance calibration data, as described in the help text immediately
% above the help for this subfunction (== help for
% 'MapLuminanceToRGBCalibrated').
%
% You would create such a file by typing, e.g.:
% btrr = 128;
% lut = the lookup table vector created by calibration...
% save 'GlobalSettings.mat' btrr lut
%
%
%
%
% Internal helper functions for Psychtoolbox - Must not be called from
% normal user code!!
%
% luttexid = PsychVideoSwitcher('GetLUTTexture', win, lut, btrr, shader);
% - Convert blue-to-luminance calibration lookup table 'lut' into a lookup
% table texture for the imaging pipeline, set it up and return a texture
% handle 'luttexid' to it. 'btrr' is the required BTRR value. 'shader' is
% the GLSL shader handle of the output formatting shader used.
%
%
% PsychVideoSwitcher(win);
% - If 'win' is a numeric onscreen window handle, perform all operations to
% implement the green channel trigger functionality for onscreen window
% 'win'. This routine uses MOGL glXXX() functions to implement drawing of
% proper trigger pixel values to the green channel for trigger creation.
%

% History:
% 05/25/08 mk  Initial incomplete implementation. Incorporates code from
%              Xiangrui Li in helper subroutines.
% 05/24/11 xl  3rd input for SwitchMode is not manditory since it is not
%              used by box version.
% 07/17/13 xl  Clarify in help text that SwitchMode only toggles between
%              two display modes.  

% GL access is needed for setup of green trigger channel in callback:
global GL;

% Cell array of trigger settings per window:
persistent triggerForWindow;

% Cell array with cached btrr and lut for window handle:
persistent cachedStuffForWindow;

% Cache info if videoswitcher is of external box type or PCI card type:
persistent VideoSwitcherIsABox;

% Subfunction dispatch:
if nargin < 1
    error('You must provide a command string!');
end

if isscalar(cmd) && isnumeric(cmd)
    % Special callback from within PTB imaging pipeline:
    if isempty(triggerForWindow) || isempty(triggerForWindow{cmd})
        % Invalid window handle or trigger disabled: Just skip.
        return;
    end

    % Extract triggerLine for this window:
    triggerSpec = triggerForWindow{cmd};
    triggerLine = triggerSpec.triggerline;

    % Trigger disabled or trigger drawn for wanted 'count' number of
    % redraws?
    if triggerLine < 0 || triggerSpec.count == 0
        % Yes. Just return.
        return;
    end
    
    % Decrement count:
    triggerSpec.count = triggerSpec.count - 1;
    triggerForWindow{cmd} = triggerSpec;
    
    % Due to a bug (or misfeature?) in imaging pipeline, we have unit 1
    % active here. Need to disable it:
    % FIXME: Fix the pipeline in next beta release cycle!
    % MK: Fixed in source code, but all Screen mex files need recompile!
    glPushAttrib(GL.ENABLE_BIT);

    glActiveTexture(GL.TEXTURE1);
    glDisable(GL.TEXTURE_RECTANGLE_ARB);
    
    % This one for unit 0 is always needed, because unit 0 is enabled by
    % design of the pipeline, not due to a bug:
    glActiveTexture(GL.TEXTURE0);
    glDisable(GL.TEXTURE_RECTANGLE_ARB);
    
    % Disable all channels for writing except green channel:
    glColorMask(GL.FALSE, GL.TRUE, GL.FALSE, GL.FALSE);
    glColor4f(0, 1, 0, 0);

    % Draw horizontal line 'triggerLine' into green channel of bound
    % framebuffer: We actually draw a stripe of 10 lines width, just to make
    % sure the VideoSwitcher picks it up, regardless of resolution.
    glBegin(GL.QUADS);
    glVertex2i(0, triggerLine);
    glVertex2i(10000, triggerLine);
    glVertex2i(10000, triggerLine+10);
    glVertex2i(0, triggerLine+10);    
    glEnd();
    
    % Reenable all channels for writing:
    glColorMask(GL.TRUE, GL.TRUE, GL.TRUE, GL.TRUE);
    glPopAttrib;
    
    % Done. Return to calling routine:
    return;
end

if ~ischar(cmd)
    error('Subcommand must be a string!');
end

if strcmpi(cmd, 'MapLuminanceToRGB')
    % Call lum2rgb() function to use simple conversion of luminance image
    % to RGB image:
    varargout{1} = lum2rgb(varargin{:});
    
    return;
end

if strcmpi(cmd, 'MapLuminanceToRGBCalibrated')
    % Call lum2rgb() function to use lut based conversion of luminance image
    % to RGB image:
    varargout{1} = lum2calrgb(varargin{:});
    
    return;
end

if strcmpi(cmd, 'SwitchMode')
    % Switch VideoSwitcher between monochrome and color mode by sending a
    % sequence of frames with special green-channel control codes:
    if nargin < 2
        error('screenId is missing!');
    end
    
    screenid = varargin{1};
    if ~ismember(screenid, Screen('Screens'))
        error('Invalid screenId provided - No such display screen!');
    end

    if nargin < 3 || isempty(varargin{2})
        enableLuminanceMode = 0; % XL: for box, this is not used
    else
        enableLuminanceMode = varargin{2};
    end
    
    if isempty(enableLuminanceMode) || ~isscalar(enableLuminanceMode) || ~isnumeric(enableLuminanceMode) || ~ismember(enableLuminanceMode, [0,1])
        error('Invalid enableLuminanceMode flag provided. Must be 0 or 1!');
    end

    if nargin < 4
        % No type parameter provided:
        if isempty(VideoSwitcherIsABox)
            % Check if type is stored in config file:
            lpath = sprintf('%sVideoSwitcherIsABox', PsychtoolboxConfigDir('VideoSwitcher'));
            if exist(lpath, 'file')
                VideoSwitcherIsABox = 1;
            else
                lpath = sprintf('%sVideoSwitcherIsACard', PsychtoolboxConfigDir('VideoSwitcher'));
                if exist(lpath, 'file')
                    VideoSwitcherIsABox = 0;
                else
                    fprintf('PsychVideoSwitcher: Could not find out if VideoSwitcher is a box or a PCI card. Just assuming it is a box...\n');
                    VideoSwitcherIsABox = 1;
                end
            end
        end
    else
        VideoSwitcherIsABox = varargin{3};
        if isempty(VideoSwitcherIsABox) || ~isscalar(VideoSwitcherIsABox) || ~isnumeric(VideoSwitcherIsABox) || ~ismember(VideoSwitcherIsABox, [0,1])
            error('Invalid VideoSwitcherIsABox flag provided. Must be 0 or 1!');
        end
    end
    
    % Enumerate all onscreen windows for screen 'screenid':
    allWindows = Screen('Windows');
    allWindowKinds = Screen('WindowKind', allWindows);
    onscreenWins = allWindows(allWindowKinds == 1);
    for win = onscreenWins
        % Is 'win' an onscreen window on screen screenid?
        if Screen('WindowScreenNumber', win) == screenid
            % This 'win'dow is an onscreen window on our target screenid.
            % Close it -- We don't want any onscreen windows on the screen
            % to switch, as they could clash/interfere with the switching
            % operation. Recycling them is impossible as we don't know
            % their configuration wrt. imaging mode, size etc.:
            Screen('Close', win);
        end
    end
    
    % Ok, our target screen 'screenid' is clear of any onscreen windows.
    % Create our own one: Fullscreen, black background clear color, no
    % imaging pipe or anything else, double-buffered. Disable all
    % sync-tests, warnings etc. for this switch action:
    oldverbosity = Screen('Preference', 'Verbosity', 1);
    oldsynclevel = Screen('Preference', 'SkipSyncTests', 2);
    oldvisuallevel = Screen('Preference', 'VisualDebugLevel', 0);
    win = Screen('OpenWindow', screenid, 0, []);
    
    % Hide mouse cursor:
    HideCursor;
    
    % Get current display size:
    [width, height]=Screen('WindowSize', win);

    % Get video refresh duration: We query nominal framerate, as this is faster
    % and sufficient for our purpose:
    ifi=Screen('NominalFramerate', win);
    if ifi == 0
        % Special case: Invalid framerate -- This usually means 60 Hz:
        ifi = 1/60;
    else
        % ifi is 1/Hz:
        ifi = 1/ifi;
    end

    if VideoSwitcherIsABox
        % Perform switching procedure for external box based VideoSwitcher:
        switchColorBox(win, width, height, ifi, enableLuminanceMode);
    else
        % Perform switching procedure for internal PCI based VideoSwitcher:
        switchColorCard(win, width, height, ifi, enableLuminanceMode);
    end
    
    % Show mouse cursor again:
    ShowCursor;

    % Close our window:
    Screen('Close', win);
    
    % Restore old settings:
    Screen('Preference', 'Verbosity', oldverbosity);
    Screen('Preference', 'SkipSyncTests', oldsynclevel);
    Screen('Preference', 'VisualDebugLevel', oldvisuallevel);
    
    return;
end

if strcmpi(cmd, 'SetTrigger')
    % Define the trigger line for a specific onscreen window:
    
    if nargin < 3
        error('You must provide the windowhandle and triggerline to the "SetTrigger" subfunction.');
    end
    
    win = varargin{1};
    if ~isscalar(win) || ~isnumeric(win) || Screen('WindowKind', win)~=1
        error('The "window" argument must be a scalar windowhandle of an onscreen window for the "SetTrigger" subfunction.');
    end

    triggerline = varargin{2};
    
    % Empty triggerline == Disable trigger:
    if isempty(triggerline)
        triggerForWindow{win} = [];
    else
        if ~isscalar(triggerline) || ~isnumeric(triggerline)
            error('The "triggerLine" argument must be a scalar numeric value for the "SetTrigger" subfunction.');
        end

        % Store this setting in our cell array:
        triggerSpec.triggerline = triggerline;
        
        % Optional count provided?
        if nargin < 4 || isempty(varargin{3})
            % Nope: Set to infinite count:
            triggerSpec.count = -1;
        else
            % Yes: Set it.
            triggerSpec.count = varargin{3};
        end
                
        triggerForWindow{win} = triggerSpec;
    end
    
    return;
end

if strcmpi(cmd, 'GetDefaultConfig')
    % Retrieve default configuration for a specific onscreen window:
    
    if nargin < 2
        screenid = 0;
    else
        % Retrieve window handle for the already opened onscreen window on
        % which VideoSwitcher is supposed to be operating:
        win = varargin{1};
        if ~isscalar(win) || ~isnumeric(win)
            error('The "window" argument must be a scalar windowhandle of an onscreen window or a scalar screenid for the "GetDefaultConfig" subfunction.');
        end

        % Screenid or windowhandle?
        if ismember(win, Screen('Screens'))
            % Screenid - Assign:
            screenid = win;
        else
            % Either windowhandle or error:
            if Screen('WindowKind', win)~=1
                error('The "window" argument must be a scalar windowhandle of an onscreen window or a scalar screenid for the "GetDefaultConfig" subfunction.');
            end

            % The 'win'dowhandle by itself is meaningless. We map the windowhandle
            % to its associated screenid and then lookup proper configuration files
            % for that screen:
            screenid = Screen('WindowScreenNumber', win);
        end
    end
    
    % First check for a per-screen file:
    lpath = sprintf('%sSettingsforScreen_%i.mat', PsychtoolboxConfigDir('VideoSwitcher'), screenid);
    if exist(lpath, 'file')
        % Per screen file: Load it.
        switcherconfig = load(lpath);
    else
        % No such per screen file. Try global default file:
        lpath = sprintf('%sGlobalSettings.mat', PsychtoolboxConfigDir('VideoSwitcher'));
        if exist(lpath, 'file')
            % Load global file:
            switcherconfig = load(lpath);
        else
            % No files at all. Setup a default config:
            switcherconfig.btrr = 128;
            switcherconfig.lut = linspace(0, 1, 257);
            fprintf('\nPsychVideoSwitcher: Warning: No meaningful configuration for VideoSwitcher available in configuration directory:\n%s\nUsing fake default settings: BTRR=128 and a linear 257 slot lut...\n\n', PsychtoolboxConfigDir('VideoSwitcher'));
            lpath = [];
        end
    end
    
    if ~isempty(lpath)
        fprintf('PsychVideoSwitcher: Info: Configuration for VideoSwitcher loaded from default config file:\n%s\n', lpath);
    end
    
    % Return BTRR setting:
    varargout{1} = switcherconfig.btrr;
    
    % Return Blue channel to luminance calibration lut, if any:
    if isfield(switcherconfig, 'lut')
        varargout{2} = switcherconfig.lut;
    else
        varargout{2} = [];
    end
    
    return;
end

if strcmpi(cmd, 'GetLUTTexture')
    if nargin < 4
        error('You must provide the window handle, lut lookup table and btrr to the "GetLUTTexture" subfunction.');
    end
    
    win = varargin{1};
    if ~isscalar(win) || ~isnumeric(win) || Screen('WindowKind', win)~=1
        error('The "window" argument must be a scalar windowhandle of an onscreen window or a scalar screenid for the "GetLUTTexture" subfunction.');
    end
    
    lut = varargin{2};
    
    if isempty(lut)
        error('You must provide the lut lookup table as 2nd argument to the "GetLUTTexture" subfunction.');
    end
    
    if ~isvector(lut)
        error('lut must be a vector!');
    end

    if size(lut, 1) ~=1
        lut = transpose(lut);
    end
    
    [channels, nslots] = size(lut);

    if channels~=1
        error('lut must be a vector!');
    end
    
    if nslots~=257
        error('lut must have 257 slots, aka columns!');
    end

    if ~isa(lut,'double')
        error('lut must by of double precision type!');
    end

    btrr = varargin{3};
    if isempty(btrr) || ~isscalar(btrr) || ~isa(btrr,'double')
        error('You must provide a valid scalar double btrr value as 3rd argument!');
    end
    
    if nargin < 5
        shader = 0;
    else
        shader = varargin{4};
        if isempty(shader)
            shader = 0;
        end
    end
    
    
    % mylut is the input lut: Row 1=Red, 2=Green, 3=Blue, 4=Alpha
    % component of a texel, 257 texels for the 257 texel wide, 1 row high
    % RGBA float texture:
    mylut = zeros(4, 257);
    
    % Fill with normalized lut:
    lut = lut / lut(257);
    
    % Hmm, a dirty trick, indeed: We shift the maximum value by a very
    % small epsilon to avoid a boundary artifact - a wrong conversion for a
    % input luminance value of exactly 1.0:
    lut(257) = lut(257) + 1e-7;
    
    % Red channel (channel 1) contains lut values, aka minimum luminance
    % values for each corresponding blue channel driver value:
    mylut(1,:) = lut;
    
    % Green channel (channel 2) contains the upper bound, ie. the channel 1
    % value of the preceding slot:
    mylut(2,1:256) = lut(2:257);

    % The upper limit of slot 257 is not defined. Set it to a value higher
    % than the minimum of slot 257:
    mylut(2,257) = lut(257);
    
    mylut(3,1:256) = repmat(btrr, 1, 256) ./ ( mylut(2,1:256) - mylut(1,1:256) );
    % This is not well-defined at all: Just set it to a high value and hope
    % for the best... As luminance values are clamped to a maximum of 1.0,
    % thi...
    mylut(3,257) = 1e10;
    
    % Build texture:
    luttex = glGenTextures(1);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, luttex);
    glTexImage2D(GL.TEXTURE_RECTANGLE_EXT, 0, GL.RGBA_FLOAT32_APPLE, 257, 1, 0, GL.RGBA, GL.FLOAT, moglsingle(mylut));

    % Make sure we use nearest neighbour sampling:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

    % And that we clamp to edge:
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
    glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);

    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);

    % Return texture handle:
    varargout{1} = luttex;

    % We store an internal backup copy of used input lut and btrr for this
    % window handle:
    wincache.btrr = btrr;
    wincache.lut = lut;
    wincache.shader = shader;
    cachedStuffForWindow{win} = wincache;
    
    return;
end

if strcmpi(cmd, 'SetBackgroundLuminanceHint')
    if nargin < 3
        error('You must provide the window handle and luminance hint value to the "SetBackgroundLuminanceHint" subfunction.');
    end
    
    win = varargin{1};
    if ~isscalar(win) || ~isnumeric(win) || Screen('WindowKind', win)~=1
        error('The "window" argument must be a scalar windowhandle of an onscreen window or a scalar screenid for the "SetBackgroundLuminanceHint" subfunction.');
    end
    
    lum = varargin{2};
    if ~isscalar(lum) || ~isnumeric(lum)
        error('The "backgroundLuminance" argument must be a scalar luminance value in range [0 ; 1] for the "SetBackgroundLuminanceHint" subfunction.');
    end
    
    % Ok, we have luminance value and window handle. Do we have relevant
    % cached info for this windowhandle?
    if isempty(cachedStuffForWindow) || isempty(cachedStuffForWindow{win})
        fprintf('PsychVideoSwitcher(''SetBackgroundLuminanceHint''): I do not have neccessary internal information for this onscreen window!\n');
        fprintf('Call PsychImaging(''OpenWindow'',...) first! This also has no function if you don''t use the calibrated VideoSwitcher driver.\n');
    end
    
    % Extract info:
    wincache = cachedStuffForWindow{win};
    
    % Only do this if the LUT based converter is used:
    if wincache.shader ~= 0
        % Map luminance value to color values via Matlab based conversion
        % routines:
        RGBPixel = lum2calrgb(lum, wincache.btrr, wincache.lut, 0);

        % Store (R,B) values in 0-1 range in red and blue channels, the
        % luminance "key" value in the green channel:
        redv = double(RGBPixel(1,1,1))/255;
        greenv = lum;
        bluev = double(RGBPixel(1,1,3))/255;

        % Find GLSL shader handle for VideoSwitcher shader:

        % Set values in shader:
        glUseProgram(wincache.shader);
        glUniform3f(glGetUniformLocation(wincache.shader, 'BackgroundPixel'), redv, greenv, bluev);
        glUseProgram(0);
    end
    
    % Done.
    return;
end

% Unknown command string - No match to any of the routines above:
error('Unknown subcommand provided!');

return; %#ok<UNRCH>

% Helper functions:
%
% Convert luminance matrix into RGB image to drive a video switcher.
%
% Usage: RGBimage = lum2rgb (lum, ratio [, trigger])
%  Inputs:
%     lum: luminance (MxN matrix with values from 0 to 1)
%     ratio: blue to red ratio of the video switcher
%     trigger: when non-zero, a trigger will be sent in current frame
%         trigger=1 or 'auto', the first line of image
%         trigger=2 or 'top',  the first line with non-zero image
%         trigger=3 or 'middle', the middle line of image
%     If you omit ratio, you should give it in this code
%  Output: RGB image (MxNx3 matrix with values from 0 to 255)
% XL, 01/03/2008
function img = lum2rgb (lum, ratio, trigger)
if nargin<3 || isempty(trigger), trigger=0; end
if nargin<2 || isempty(ratio), ratio=128; end
if nargin<1, disp('Usage: RGBimage = lum2rgb (lum, ratio [, trigger])'); return; end

sz=size(lum);
img = uint8(zeros(sz(1),sz(2),3));   % preallocate memory
lum=lum*255;	% now 0 to 255
img(:,:,3) = uint8(min(255,floor((ratio+1)/ratio*lum)));       % first calculates BLUE
img(:,:,1) = uint8((ratio+1)*lum-double(img(:,:,3))*ratio);	% then remainder is achieved by RED

if trigger==0, return; end

switch trigger
    case {2, 'top'}
        line = 1;    % first line of the image
    case {3, 'middle'}
        line = round(sz(1)/2);    % middle of image
    otherwise
        line = ceil(find(lum',1,'first')/sz(2));    % the first line where lum is non-zero
end
img(line,:,2) = uint8(255);         % set green to 255. Will send trigger when displayed
return;

% Convert luminance matrix into RGB image to drive a video switcher.
% This uses the calibrated luminance table saved in calibratedlum.mat.
% For detail, check the paper http://lobes.usc.edu/Journals/JNM03.pdf
% 
% The steps to get the table:
%  1. Switch the video switcher to grayscale mode;
%  2. Set up the equipement to accurately measure screen luminance;
%     you can use a photometer or data acquisition system;
%  3. Measure 257 luminance levels at RGB of [0 0 b] and [btrr 0 255],
%     where b is 0:255, and btrr the blue to red ratio of switcher; 
%  4. Store 257 luminance in variable callum and normalize them:
%     callum=callum/callum(257);
%  5. Save it under a proper folder with Matlab path: 
%     save calibratedlum callum;
% 
% Usage: RGBimage = lum2calrgb (lum, ratio [, trigger])
%  Inputs:
%     lum: luminance (MxN matrix with values from 0 to 1)
%     ratio: blue to red ratio of the video switcher
%     trigger: when non-zero, a trigger will be sent in current frame
%         trigger=1 or 'auto', the first line of image
%         trigger=2 or 'top',  the first line with non-zero image
%         trigger=3 or 'middle', the middle line of image
%     If you omit ratio, you should give it in this code
%  Output: RGB image (MxNx3 matrix with values from 0 to 255)
% 
% Note: this runs slow when lum matrix is large. 200x200 can take 1 second.
% XL, 01/08/2008
function img = lum2calrgb (lum, ratio, callum, trigger)
if nargin<4 || isempty(trigger), trigger=0; end
if nargin<3, disp('Usage: RGBimage = lum2calrgb (lum, ratio, lut [, trigger])'); return; end

if length(callum)~=257
    error('Calibrated luminance table must have 257 entries.');
end

sz=size(lum);
blue=ones(sz)*255;
for i=1:256
    ind=find(lum>=callum(i) & lum<callum(i+1));
    blue(ind)=i-1; %#ok<FNDSB>
end
img=zeros(sz(1),sz(2),3);   % preallocate memory
img(:,:,3)=blue;
img(:,:,1)=(lum-callum(blue+1))./(callum(blue+2)-callum(blue+1))*ratio;  % red
img=uint8(img);

if trigger==0, return; end
switch trigger
    case {2, 'top'}
        line = 1;    % first line of the image
    case {3, 'middle'}
        line = round(sz(1)/2);    % middle of image
    otherwise
        line = ceil(find(lum',1,'first')/sz(2));    % the first line where lum is non-zero
end
img(line,:,2) = uint8(255);         % set green to 255. Will send trigger when displayed

return;

% Switch between normal color mode and high grayscale mode of a box video switcher. 
% Usage for box switcher:
%  switchColor(windowOrScreenNumber) to switch back and forth between two modes.
% Call it to switch to grayscale mode after you open a window, and 
% call it again to switch back to normal color mode after you finish display.
% XL, 01/08/2008
function switchColorBox(w, width, height, ifi, enableLuminanceMode) %#ok<INUSD>

% Doesn't help doesn't hurt either:
Screen('FillRect',w,0);

ppb=1/ifi*width*height/4e6; % pixels per bit based on 4 MHz detection frequency

% these 2 signals will make a switch
bits=[1 1 1 0 0 1 0 1 0 1 0 1 1 1 1 0 0 0 1 1 0 1 1 0 1; ... 
      1 1 0 1 1 0 1 1 1 0 1 1 0 1 1 0 1 0 1 0 1 1 0 0 1]; 

% factor means the percent of time to display pixels within a frame.
% Here we use a large range of factor so it will work for different settings.
range=[0.6 0.85];   % you can narrow the range at specific monitor setting
step=0.2/ppb; % in case of switch failure, reduce 0.2 to about 0.15
pixelPos=round((0:25)*ppb/range(1));
img1=zeros(1, pixelPos(26), 3); % one line image, will change green layer
img2=img1;
for i=1:25
    img1(1, pixelPos(i)+1:pixelPos(i+1), 2)=bits(1,i)*255;
    img2(1, pixelPos(i)+1:pixelPos(i+1), 2)=bits(2,i)*255;
end
tex(1)=Screen('MakeTexture',w,uint8(img1));
tex(2)=Screen('MakeTexture',w,uint8(img2));

for i=1:2
    for factor=(range(1):step:range(2))/range(1)
        rd=rand*ppb*0.5; % random shift up to half a bit
        Screen('DrawTexture',w, tex(i),[],[rd 0 pixelPos(26)/factor+rd 1]); 
        Screen('Flip',w);
    end
end

% Release textures:
Screen('Close', tex);

return;

% Switch between normal color mode and high grayscale mode of a card video switcher. 
% Usage: switchColorCard [(color,windowOrScreenNumber)]
%   Optional inputs: 
%     color: 0, 'BW' or 'gray' to switch to high grayscale monochrome mode, 
%        otherwise (or omit) to switch to color mode.
%     windowOrScreenNumber: windowPtr or screen number. Default is current onscreen.
% XL, 11/2007

function switchColorCard(w, width, height, ifi, enableLuminanceMode)

Screen('FillRect',w,0);

ppb=1/ifi*width*height/4e6; % pixels per bit based on 4 MHz detection frequency

% 1st row to gray, 2nd color
bits=[1 1 1 0 0 1 0 1 0 1 0 1 1 1 1 0 0 0 1 1 0 1 1 0 1; ... 
      1 1 0 1 1 0 1 1 1 0 1 1 0 1 1 0 1 0 1 0 1 1 0 0 1]; 

if enableLuminanceMode==1
    % Switch to monochrome mode:
    color=1;
else
    % Switch to RGB mode:
    color=2;
end

% factor means the percent of time to display pixels within a frame.
% Here we use a large range of factor so it will work for different settings.
range=[0.6 0.85];   % you can narrow the range at specific monitor setting
step=0.2/ppb; % in case of switch failure, reduce 0.2 to about 0.15
pixelPos=round((0:25)*ppb/range(1));
img=zeros(1, pixelPos(26), 3); % one line image, will change green layer
for i=1:25
    img(1, pixelPos(i)+1:pixelPos(i+1), 2)=bits(color,i)*255;
end
tex=Screen('MakeTexture',w,uint8(img));

for factor=(range(1):step:range(2))/range(1)
    rd=rand*ppb*0.5; % random shift up to half a bit
    Screen('DrawTexture',w,tex,[],[rd 0 rd+pixelPos(26)/factor 1]);
    Screen('Flip',w); 
end

% MK: Release textures:
Screen('Close', tex);

return;
