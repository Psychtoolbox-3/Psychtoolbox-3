function varargout = PsychDataPixx(cmd, varargin)
% PsychDataPixx - High level control driver for the VPixx - DataPixx device.
%
% This driver provides common high-level functionality for interaction with
% the VPixx Technologies DataPixx device. The driver provides high-level
% functions for basic device operations and for timestamping of visual
% stimulus onset. It also provides functionality needed for the
% PsychImaging() command and the Screen() image processing pipeline to
% properly setup the high-precision video display modes and stereo display
% modes of the DataPixx.
%
% The driver is dependend on the low-level MEX file Datapixx() driver to be
% present. It will fail if that MEX file is missing or dysfunctional.
%
% Normally you will have a DataPixx connected when using this driver, but
% for basic code development and testing you can call
% PsychDataPixx('SetDummyMode', 1); and then run without a connected
% DataPixx device.
%
% For control of visual stimulus display on the DataPixx, read the relevant
% sections of the help of PsychImaging. That command will take care of most
% of the details of visual setup.
%
%
% Commands and their syntax:
% --------------------------
%
% oldmode = PsychDataPixx('SetDummyMode' [, mode]);
% - Switch driver into dummy mode if 'mode' is set to 1. In dummy mode, the
% code mostly operates without a DataPixx device connected. This is useful
% for basic debugging of code without access to the device. The simulated
% device returns more or less meaningful values, good enough for initial
% code development.
%
% A 'mode' setting of zero enabled normal mode of operation.
% You must call this command before any other command for proper operation!
%
%
% PsychDataPixx('Open');
% - Open a connection to the device, with default settings. If the device
% is already open then this will do almost nothing, otherwise the
% connection is opened.
%
% You must call this function before all other remaining commands!
%
%
% PsychDataPixx('Close');
% - Close a connection to the device. Multiple virtual connections can be
% simultaneously open. If this function call closes the last virtual
% connection, then the real physical connection to the device is closed and
% the driver is reset.
%
% If you used PsychImaging() to make use of Datapixx graphics mode, it will
% automatically close that connection once the corresponding onscreen
% window is closed. As a simple rule, you only need to call the 'Close'
% function if you also called the 'Open' function before.
%
% Call this command as last command in your script! All further
% PsychDataPixx() commands will be invalid after a 'Close' call!
%
%
% status = PsychDataPixx('GetStatus' [,newstatus]);
% - Retrieve a struct 'status' with the complete driver internal state.
% Optionally assign new state 'newstatus' to driver. Only assign a new
% state if you *really* know what you're doing!
%
%
% oldmode = PsychDataPixx('LogOnsetTimestamps', mode);
% - Return current timestamp acquisition mode as optional return argument
% 'oldmode'. Set a new mode, according to argument 'mode'. 'mode' can
% be one of:
%
% 0 = Disable all timestamp acquisition and logging.
% 1 = Log a visual stimulus onset timestamp for the next Screen('Flip') or
%     Screen('AsyncFlipBegin') command after this function call, then stop
%     logging again. This is a one-shot timestamp function.
% 2 = Log all Screen('Flip') and Screen('AsyncFlipBegin') visual stimulus
%     onset timestamps until this function is called again with a 'mode'
%     setting of 1 or 0.
%
% Timestamping is disabled by default (mode == 0), as it incurs a bit of
% computational overhead to acquire and log timestamps, typically up to 2-3
% msecs of extra time per 'Flip' command.
%
% The driver will log all timestamps in a buffer, which can be read out via
% the following timestamp related commands.
%
%
% [boxTime, getsecsTime, flipCount, currentFlipcount] = PsychDataPixx('GetLastOnsetTimestamp');
% - Retrieve visual stimulus onset timestamps related to last logged
% "Flip". This call will only return valid values after a Flip command is
% finished, ie., after a successfull call to Screen('Flip') or
% Screen('AsyncFlipCheckEnd') or Screen('AsyncFlipEnd'), and only if
% logging of timestamps was enabled by setting the 'LogOnsetTimestamps'
% mode (see above) to a non-zero value. It always returns the most recently
% acquired timestamp.
%
% 'boxTime' is the time in seconds of stimulus onset as measured by the
% DataPixx internal clock.
% 'getsecsTime' is the time in seconds of stimulus onset as measured by the
% Psychtoolbox GetSecs() clock.
% 'flipCount' is the "serial number" of the Screen flip command to which
% the returned timestamps refer.
% 'currentFlipcount' is the "serial number" of the most recent Screen flip
% command.
%
% Computation of 'getsecsTime' will become inaccurate if you don't
% frequently call PsychDataPixx('GetPreciseTime') to resynchronize the
% clocks! A more accurate approach to retrieve GetSecs() style timestamps
% is to use the PsychDataPixx('GetTimestampLog', 1) function to retrieve
% all logged timestamps at the end of a session and apply some
% high-precision remapping of timestamps.
%
%
% [getsecs, boxsecs, confidence] = PsychDataPixx('GetPreciseTime');
% - Query both the DataPixx clock and GetSecs clock and compute which GetSecs
% time 'getsecs' corresponds to which Datapixx time 'boxsecs' and how
% reliable this correspondence is in 'confidence'. 'confidence' is the
% margin of error between both timestamps.
% This function implicitely synchronizes both clocks to get more precise
% values for 'getsecstime' from calls to PsychDataPixx('GetLastOnsetTimestamp');
% The function is automatically called once at device open time.
%
%
% log = PsychDataPixx('GetTimestampLog' [, remapPrecise = 0]);
% - Fetch the full log of all acquired visual stimulus onset timestamps
% into the variable 'log'. 'log' is a 3-by-n matrix for the n timestamps,
% one column encodes one stimulus onset, aka one logged invocation of
% Screen('Flip') et al.
%
% log(1,i) = i'th timestamp as measured via the Datapixx internal clock.
% log(2,i) = i'th timestamp as measured via the Psychtoolbox GetSecs clock.
% log(3,i) = Flip serial number of i'th sample.
%
% If the optional parameter 'remapPrecise' is set to 1, then an expensive
% calibration procedure is used to compute very precise GetSecs()
% timestamps for row 2 of the returned matrix, otherwise a cheap but
% inaccurate remapping is performed. The expensive procedure can take up to
% 1 second, but is highly recommended for precise timestamps. See
% explanation of 'BoxsecsToGetsecs' function below for more details. You'll
% typically call this function only at the end of an experiment session.
%
%
% PsychDataPixx('ClearTimestampLog');
% - Clear the log with all currently acquired timestamps.
%
%
% tgetsecs = PsychDataPixx('FastBoxsecsToGetsecs', t);
% - Map given timestamp 't' in Datapixx clock time to Psychtoolbox GetSecs
% time and return it in 'tgetsecs'. This mapping is only precise if you
% call PsychDataPixx('GetPreciseTime'); frequently. For a more accurate
% remapping use the following PsychDataPixx('BoxsecsToGetsecs', t);
% instead, which is more time consuming to execute though.
%
%
% [tgetsecs, sd, ratio] = PsychDataPixx('BoxsecsToGetsecs', t);
% - Perform remapping of a vector of Datapixx timestamps 't', ie.,
% timestamps expressed in Datapixx clock time, into Psychtoolbox GetSecs
% timestamps 'tgetsecs'. Return measures of accuracy of remapping.
%
% This must be called while the device is still open! It performs the same
% calibration procedure as PsychDataPixx('GetTimestampLog', 1); and can
% easily take up to 1 second of time.
%
% The remapping is performed by fitting a mapping function to the set of
% acquired timestamp samples from both the Psychtoolbox GetSecs clock and
% the DataPixx clock. The best fit function is used for mapping DataPixx
% timestamps in 't' to PTB timestamps in 'tgetsecs'. 'sd' is the standard
% deviation of the best-fit function wrt. to fitting error. 'ratio' is an
% estimate of the clock speed difference between the host clock and the
% Datapixx clock.
%
%
% PsychDataPixx('RequestPsyncedUpdate');
% - Request an update of the DataPixx register block in sync with the next
% visual stimulus onset as triggered by Screen('Flip') or
% Screen('AsyncFlipBegin'). This will emit the neccessary PSYNC pixel
% sequence at next flip and tell the device to apply all pending settings
% on reception of that PSYNC token.
%
%
% count = PsychDataPixx('FlipCount');
% - Return "serial number" of the last executed visual stimulus update via
% one of the Screen "Flip" commands.
%
%
% PsychDataPixx('ExecuteAtFlipCount', targetFlipCount, commandString);
% - Request execution of 'commandString' via eval() function in sync with
% the Screen('Flip') or Screen('AsyncFlipBegin') command which will
% cause the stimulus onset with "serial number" 'targetFlipCount'. You can
% ask for the current "serial number" by calling current = PsychDataPixx('FlipCount');
% Providing an empty 'targetFlipCount' means to execute at next flip.
%
% Datapixx related commands in 'commandString' will be scheduled for
% execution on the device at next visual stimulus onset by use of the psync
% mechanism.
% 
% Not all commands are allowed inside 'commandString'! Any call to Screen()
% or any OpenGL command is forbidden. Datapixx commands which could block
% execution to wait for the Datapixx , e.g., Datapixx('RegWrRd') are
% problematic, and the calls Datapixx('RegWrPixelSync') and
% Datapixx('RegWrRdPixelSync') are strictly forbidden!
%
% This function is mostly used by other higher-level DataPixxToolbox
% functions, e.g., for I/O (sound, digital or analog), in order to
% synchronize their operation and schedules to visual stimulus updates.
%
%
% PsychDataPixx('LoadIdentityClut', win);
% - Load an identity CLUT into the device at next stimulus onset, ie., next
% 'Flip' command for window 'win'. This is a pure convenience function for
% that common case. Normally you'd use Screen('LoadNormalizedGammaTable', win, clut, 2);
% to upload a new 'clut' table at next 'Flip' command execution.
%
%
% oldverbosity = PsychDataPixx('Verbosity' [,verbosity]);
% - Retrieve and optionally set a new level of 'verbosity' for driver debug
% output. verbosity can be 0 for no output, 1 for only errors, 2 for
% additionally warnings, 3 for additional info, 4 for very verbose output.
%
%
% oldtimeout = PsychDataPixx('PsyncTimeoutFrames' [, timeoutFrames]);
% -- Query and/or set timeout (in video refresh cycles) for recognition of
% PSYNC codes by the device. By default, the timeout is set to the
% equivalent of 5 minutes, which means: If the device is instructed to wait
% for a PSYNC marker token in the video stream, but doesn't receive the
% expected token within 5 minutes, it will continue processing as if the
% token was received. This is used to prevent device-hangs on programming
% errors or other malfunctions. The default setting of 5 minutes is very
% generous. You can override this default and set an arbitrary timeout
% between 1 and 65535 video refresh cycles with this function.
%
%
%
%
% Normally you won't call the following functions yourself, as Psychtoolbox
% automatically performs neccessary setup during calls to PsychImaging()
% that are related to the DataPixx device ("help PsychImaging" for
% details).
%
%
% PsychDataPixx('SetVideoMode', mode);
% -- Switch DPixx device immediately into video processing mode 'mode'.
% Type "Datapixx SetVideoMode?" for a list of available 'mode' settings and
% explanation of their purpose and behaviour.
%
%
% PsychDataPixx('SetVideoHorizontalSplit', mode);
% -- Switch DPixx device immediately into video split mode 'mode'.
% Type "Datapixx SetVideoHorizontalSplit?" for a list of available 'mode'
% settings and explanation of their purpose and behaviour.
%
%
% PsychDataPixx('SetVideoVerticalStereo', mode);
% -- Switch DPixx device immediately into vertical stereo video mode 'mode'.
% Type "Datapixx SetVideoVerticalStereo?" for a list of available 'mode'
% settings and explanation of their purpose and behaviour.
%
%
% PsychDataPixx('EnableVideoScanningBacklight');
% -- Enable ViewPixx panels scanning backlight.
%
%
% PsychDataPixx('DisableVideoScanningBacklight');
% -- Disable ViewPixx panels scanning backlight.
%
%
% PsychDataPixx('EnableVideoStereoBlueline');
% -- Enable detection and handling of blue-line stereo sync lines by the
% DataPixx.
%
%
% PsychDataPixx('DisableVideoStereoBlueline');
% -- Disable detection and handling of blue-line stereo sync lines by the
% DataPixx.
%
%
% Internal commands: NOT FOR USE BY PURE MORTALS!
%
% PsychDataPixx('PerformPostWindowOpenSetup');
% -- Called by PsychImaging('OpenWindow') after opening a window on the
% device. Performs all low-level setup for use of DataPixx.
%
%
% rc = PsychDataPixx('CheckGPUSanity', window, xoffset);
% -- Perform online-test of GPU identity gamma tables and DVI-D display
% encoders. Try to correct problems with wrong identity gamma tables and at
% least detect problems with (spatio-)temporal display dithering. Returns
% rc == 0 on full success, rc > 0 on failure.
%
%
% PsychDataPixx(-1);
% -- PTB callback: Request immediate emission of a PSYNC'ed RegWrRd command
% to driver if a PSYNC is pending. Otherwise noop.
%
%
% PsychDataPixx(0);
% -- PTB callback at 'Flip' time: Do whatever DPixx related work is
% pending for next 'Flip', e.g., DPixx clut updates, writing of PSYNC pixel
% sequences, preparation of logging of timestamps etc.
%
%
% PsychDataPixx(1, clut);
% -- Callback for PTB imaging pipeline: Ultrafast upload of given 'clut'
% into the device at next bufferswap.
%
%
% PsychDataPixx(2);
% -- Callback for PTB imaging pipeline after successfull finish of a
% bufferswap. Do whatever needs to be done, e.g., incrementing the
% swapcounter and logging of timestamps.
%
%
% PsychDataPixx(3);
% -- Callback at time of closing the DataPixx onscreen window: Perform
% whatever actions are needed at shutdown, e.g., restoring identity clut on
% device, restoring default video mode and settings, closing device
% connection, cleanup of data structures etc.
%
%

% Possible improvements:
%
% * More options for different subfunctions? More subfunctions for
%   high-level stuff? Response box functionality should probably go into
%   another high level driver which can simply utilize the timestamp
%   remapping routines inside this driver?

% History:
% 19.12.2009  mk  Written. First iteration.
% 20.01.2010  mk  Switch psync-command to preflip operations: Now done as
%                 part of preflip ops, potentially long before the swap
%                 deadline, not from the prebufferswap callback anymore.
%                 This should improve timing precision of scheduled
%                 stimulus onset and allows to use PSYNC mechanism also for
%                 async-flips, because now it is thread-safe. We now make
%                 use of the 16-bit timeout value for Psync and set it to a
%                 default 5 minutes, but possible to override by usercode.
%                 We also setup the device to only accept psync in the
%                 first scanline of the display and to blank that scanline
%                 to black.
%                 --> No remaining limitations for use with PTB and no
%                 visual artifacts from psync mechanism :-)
%
% 07.09.2010  mk  Implement online diagnostic and correction for proper
%                 function of GPU identity gamma tables and encoders.
%                 Called as 'CheckGPUSanity' from imaging pipeline setup
%                 code in BitsPlusPlus() for high precision display mode.
%
% 27.11.2011  mk  Improve online correction routine for LUT's. Need to
%                 handle GPUs with more than 256 slot LUT's sensibly. E.g.,
%                 the NVidia QuadroFX-3800 has 2048 LUT slots.
%
% 11.09.2012  mk  Add support for stereo sync line handling, and for
%                 scanning backlight control of the ViewPixx.

% Need GL constant for low-level OpenGL calls. Already initialized by
% PsychImaging() at first time of invocation of the driver:
global GL;

% Global dpx device data structure with all our variables:
global dpx;

% Cold init at full restart:
if isempty(dpx)    
    dpx.verbosity = 1;
    dpx.refcount = 0;
    dpx.needpsync = 0;
    dpx.needonsettimestamp = 0;
    dpx.dummymode = 0;
    dpx.lasttOnset = nan;
    dpx.lastLoggedOnsetFramecount = 0;
    dpx.swapcount = 0;
    dpx.logonsettimestamps = 0;
    dpx.timestampLogCount = 0;
    dpx.onsetTimestampLog = [];
    dpx.timeout = [];

    % Allocate 2 hours worth of timestamps at 200 Hz flip rate:
    dpx.timestampLogPreallocSize = 2 * 3600 * 200;

    % Syncmethod 1 is empirically the most accurate for clock sync:
    dpx.syncmode = 1;
    % Use a maximum duration for syncClocks call of 0.5 seconds.
    dpx.maxDuration = 0.5;

    % No optimum for abort: Try the best you can!
    dpx.optMinwinThreshold = 0.0;
    
    % Maximum allowable minwin is 1.3 msecs, the expected maximum for an ok
    % but not great USB write transaction:
    dpx.maxMinwinThreshold = 0.0013;
    
    % Logging array for host <-> device clock sync samples:
    dpx.syncSamples = [];
    
    % Clockratio is assumed to be 1, until we decide to implement
    % clockRatio calibration routines:
    dpx.clkRatio = 1.0;
    
    % ID unknown until 'Open' call:
    dpx.ID = '';
    
    % Datapixx already open?
    if doDatapixx('IsReady')
        % Yup: Some external client holds an open connection. Artificially
        % increase our refcount to prevent us from accidentally closing
        % that connection down when we close our connections:
        dpx.refcount = 1;
    end
    
    % Init flip-schedule to empty, ie., one free slot. It will
    % automatically grow if more than 1 slot is ever needed:
    dpx.targetswapcounts = -1;
    dpx.scheduledcommands = cell(1,1);
    
    % No windowhandle:
    dpx.window = [];    
end

if cmd == 0
    % Fast callback from within Screen - This one is issued each time
    % before scheduling a doublebuffer-swap on the Datapixx display device.
    % Done as part of preflip operations on the main thread. There can be
    % an arbitrarily long unknown delay between invocation of this callback
    % and the actual bufferswap, but there will be no graphics related
    % activity at all for the display between return from this callback and actual
    % bufferswap -- the backbuffer is in its final state, ready for display
    % at next bufferswap.
    %
    % The callback has to perform any work scheduled for the next stimulus
    % onset that is related to the Datapixx:

    if isempty(GL)
        fprintf('PsychDataPixx: FATAL ERROR! OpenGL struct GL not initialized in Screen callback!! Driver bug?!?\n');
        return;
    end

    % Find all work-to-be-done items for target swapcount of next flip:
    worklets = find(dpx.targetswapcounts == dpx.swapcount + 1);
    if ~isempty(worklets)
        % Iterate over all work items on the todo-list:
        for workletidx = worklets
            % Execute!
            eval(dpx.scheduledcommands{workletidx});
            
            % Done. Free slot in schedule:
            dpx.scheduledcommands{workletidx} = [];
            dpx.targetswapcounts(workletidx) = -1;
        end
        
        % Need PSYNC'ed submission of these commands to DataPixx:
        dpx.needpsync = 1;
    end
    
    % We always unconditionally draw scanline 1 in all-black, so PSYNC may
    % work reliably:
    glRasterPos2i(0, 1);
    glDrawPixels(size(dpx.blackline, 2), 1, GL.RGB, GL.UNSIGNED_BYTE, dpx.blackline);
    
    % Any work to do in sync with next stimulus onset?
    if dpx.needpsync || dpx.needonsettimestamp
        % Yes! Previous calls to PsychDataPixx have emitted Datapixx()
        % register status update commands. These need to be applied on
        % reception of the next psync token as part of stimulus onset.
        
        if dpx.needonsettimestamp
            % Emit latch command for timestamping:
            doDatapixx('SetMarker');
        end
        
        % Define psync code:
        dpx.psync = [ 255 0 0 ; 0 255 0 ; 0 0 255 ; 255 255 0 ; 255 0 255 ; 0 255 255 ; 255 255 255 ; 0 0 0 ]'; 

        % Make code unique for this frame, so no two consecutive psync
        % updates ever have the same psync token:
        % bitand(bitshift(x, -0), 255)
        dpx.psync(2, 8) = mod(dpx.swapcount + 1, 256);
        
        % Blit it to top-left corner of framebuffer:
        glRasterPos2i(10, 1);
        glDrawPixels(size(dpx.psync, 2), 1, GL.RGB, GL.UNSIGNED_BYTE, uint8(dpx.psync));

        % Following is asyncflip-safe, as opposed to the old implementation
        % which would only work for syncflips:

        % Emit psync - synced register writecommand:
        % Set device timeout for PSYNC detection to dpx.timeout video
        % refresh cycles:
        doDatapixx('RegWrPixelSync', dpx.psync, dpx.timeout);

        % Set "psync work to do" flag to "psync pending"
        dpx.needpsync = 2;
    end
    
    % Return control to Screen():
    return;
end

if cmd == 1
    % Fast callback from within Screen's imaging pipeline for setting a new
    % hardware clut in device at next doublebufferswap:
    
    % Upload clut to device:
    doDatapixx('SetVideoClut', varargin{1});
    
    % Request PSYNC'ed application at next stimulus onset:
    dpx.needpsync = 1;

    if dpx.verbosity > 5
        fprintf('CLUT UPDATE!\n\n');
        disp(varargin{1});
        fprintf('------------\n\n');
    end

    % Return control to Screen():
    return;
end

if cmd == 2
    % Fast callback from within Screen's imaging pipeline some arbitrary
    % time after successfull completion of a bufferswap, but certainly
    % before any other graphics related activity for the onscreen window
    % will happen. This is called at end of execution of Screen('Flip') or
    % Screen('Asyncflipend') or successfull completion of
    % Screen('Asyncflipcheckend'):

    % Reset "psync work to do" flag to "none":
    dpx.needpsync = 0;

    % Increment swap-completed count:
    % Could do this to query for true flip count: dpx.swapcount = varargin{1}
    % But this is better for now, as the above doesn't start with value 1
    % at our first invocation, but is already at around 7 due to prior
    % flips of setup code - Avoid the hassle, keep it simple:
    dpx.swapcount = dpx.swapcount + 1;
    
    % Were there any onset timestamps requested for the finished swap?
    if dpx.needonsettimestamp
        % Yes. Get'em from device:
        doDatapixx('RegWrRd');
        dpx.lasttOnset = doDatapixx('GetMarker');
        dpx.lastLoggedOnsetFramecount = dpx.swapcount;
        
        % Was this a singleshot timestamp request (==1) or a continuous one
        % (==2)
        if dpx.needonsettimestamp == 1
            % Single-Shot: Reset to "don't need timestamp":
            dpx.needonsettimestamp = 0;
        end
        
        % Timestamp logging requested?
        if dpx.logonsettimestamps
            dpx.timestampLogCount = dpx.timestampLogCount + 1;
            dpx.onsetTimestampLog(:, dpx.timestampLogCount) = [dpx.lasttOnset ; box2GetSecsTime(dpx.lasttOnset); dpx.swapcount];
        end
    end
    
    % Return control to Screen():
    return;
end

if cmd == 3
    % Callback from within Screen at window close time as part of
    % Screen('Close', window), Screen('CloseAll') or auto-close on error or
    % Screen shutdown. Called with the OpenGL subsystem already offline:
    
    % Just call the 'ResetOnWindowClose' routine for now:
    cmd = 'ResetOnWindowClose';
end

if strcmpi(cmd, 'ResetOnWindowClose')
    % Shutdown and cleanup operations when the onscreen window is closed
    % that displays on the Datapixx:
    if dpx.verbosity > 3
        fprintf('PsychDataPixx: Closing device connection for this onscreen window...\n');
    end
    
    % Immediately load a standard identity CLUT into the device:
    linear_lut =  repmat(linspace(0, 1, 256)', 1, 3);
    doDatapixx('SetVideoClut', linear_lut);
    
    % Reset videomode to pass-through:
    doDatapixx('SetVideoMode', 0);
    
    % Set stereomode to "auto" and horizontal split to "auto":
    doDatapixx('SetVideoHorizontalSplit', 2);
    doDatapixx('SetVideoVerticalStereo', 2);
    
    % Disable Pixelsyncline:
    doDatapixx('SetVideoPixelSyncLine', 0, 0, 0);
    
    % Apply all changes immediately:
    doDatapixx('RegWrRd');

    % Delete windowhandle:
    dpx.window = [];
    
    % Perform standard close op if needed:
    cmd = 'Close';
end

if strcmpi(cmd, 'Verbosity')
    % Return old verbosity:
    varargout{1} = dpx.verbosity;
    
    if nargin >= 2 && ~isempty(varargin{1})
        % Assign new verbosity:
        dpx.verbosity = varargin{1}; %#ok<NASGU>
    end
    
    return;
end

if strcmpi(cmd, 'PsyncTimeoutFrames')
    % Return old timeout in frames:
    varargout{1} = dpx.timeout;
    
    if nargin >= 2 && ~isempty(varargin{1})
        % Assign new verbosity:
        dpx.timeout = varargin{1}; %#ok<NASGU>
    end
    
    return;
end

if strcmpi(cmd, 'LogOnsetTimestamps')
    % Switch state of automatic timestamp logging:
    
    % Return old setting:
    varargout{1} = dpx.logonsettimestamps;
    
    if nargin < 2 || isempty(varargin{1})
        error('PsychDataPixx: LogOnsetTimestamps: You must provide the "enable" flag to subfunction!');
    end
    
    % Assign new mode:
    dpx.logonsettimestamps = varargin{1}; %#ok<NASGU>
    
    % Set new prealloc size for timestamp log?
    if nargin >= 3 && ~isempty(varargin{2})
        % Yes: Assign it.
        dpx.timestampLogPreallocSize = varargin{2};

        % Clear current log, so it gets realloced below:
        dpx.onsetTimestampLog = [];
        dpx.timestampLogCount = 0;
    end
    
    % Disable logging requested?
    if dpx.logonsettimestamps == 0
        % Yes. Reset the need timestamp flag to '0' aka off:
        dpx.needonsettimestamp = 0;
    end

    % Singleshot logging requested, ie., only for next Screen('Flip')?
    if dpx.logonsettimestamps == 1
        % Yes. Set the need timestamp flag to '1' aka singleshot:
        dpx.needonsettimestamp = 1;
    end
    
    % Continouous logging requested, ie., for each Screen('Flip')?
    if dpx.logonsettimestamps == 2
        % Yes. Set the need timestamp flag to '2' aka always:
        dpx.needonsettimestamp = 2;
    end
    
    % Allocate timestamp log if not already alloated:
    if dpx.logonsettimestamps && isempty(dpx.onsetTimestampLog)
        % Always allocate dpx.timestampLogPreallocSize number of elements:
        dpx.onsetTimestampLog = zeros(3, dpx.timestampLogPreallocSize);
    end

    return;
end

if strcmpi(cmd, 'GetPreciseTime')
    % Perform a clocksync between host and device, return triplet of
    % hosttime, corresponding DataPixx device clock time and a confidence
    % interval for the accuracy of clock sync. Return all three:
    res = syncClocks;
    varargout{1} = res(1);
    varargout{2} = res(2);
    varargout{3} = res(3);
    return;
end

if strcmpi(cmd, 'GetStatus')
    % Return copy of internal driver status structure:
    varargout{1} = dpx;
    
    % Optionally assign driver status from external input argument:
    if nargin > 1 && ~isempty(varargin{1})
        dpx = varargin{1};
    end
    
    return;
end

if strcmpi(cmd, 'ClearTimestampLog')
    % Clear/Delete the current timestamp log of flip onset timestamps:
    dpx.onsetTimestampLog = [];
    dpx.timestampLogCount = 0;
    return;
end

if strcmpi(cmd, 'GetTimestampLog')
    if nargin >= 2 && ~isempty(varargin{1})
        if varargin{1} > 0
            % Precision timestamps in GetSecs() time requested:
            % We perform a - expensive - high precision remapping of the
            % raw timestamps in row 1 of the current log and store the
            % resulting "best fit" GetSecs timestamps in row 2:
            [timing, sd, clockratio] = box2GetSecsTimePostHoc(dpx.onsetTimestampLog(1, 1:dpx.timestampLogCount));
            dpx.onsetTimestampLog(2, 1:dpx.timestampLogCount) = timing;
            varargout{2} = sd;
            varargout{3} = clockratio;
        end
    end
    
    % Get a copy of the current timestamp log of flip onset timestamps:
    varargout{1} = dpx.onsetTimestampLog(:, 1:dpx.timestampLogCount);
    return;
end

if strcmpi(cmd, 'FastBoxsecsToGetsecs')
    if nargin < 2 || isempty(varargin{1})
        error('PsychDataPixx: FastBoxsecsToGetsecs: Vector with timestamps to remap missing!');
    end

    % Return result of fast remapping, based on last syncClocks run:
    varargout{1} = box2GetSecsTime(varargin{1});
    return;
end

if strcmpi(cmd, 'BoxsecsToGetsecs')
    if nargin < 2 || isempty(varargin{1})
        error('PsychDataPixx: BoxsecsToGetsecs: Vector with timestamps to remap missing!');
    end

    [timing, sd, clockratio] = box2GetSecsTimePostHoc(varargin{1});
    varargout{1} = timing;
    varargout{2} = sd;
    varargout{3} = clockratio;

    return;
end

if strcmpi(cmd, 'GetLastOnsetTimestamp')
    % Get a copy of the onset flip timestamp of the most recent flip:
    varargout{1} = dpx.lasttOnset;
    
    % Also provide time remapped to GetSecs time, based on last
    % syncClocks() run:
    varargout{2} = box2GetSecsTime(dpx.lasttOnset);

    % Also return the specific swapcount to which the timestamp refers:
    varargout{3} = dpx.lastLoggedOnsetFramecount;

    % Also return the current swapcount of last flip:
    varargout{4} = dpx.swapcount;
    return;
end

if strcmpi(cmd, 'RequestPsyncedUpdate')
    % External calling code needs a psync - synced register update at next
    % Screen('Flip')'ed stimulus onset:
    if dpx.needpsync > 1
        % Ouch! Request happens after PSYNC command emission. This means
        % chaos! We error out:
        fprintf('PsychDataPixx: ERROR: RequestPsyncedUpdate called inbetween a call to\n');
        fprintf('PsychDataPixx: ERROR: Screen("Drawingfinished/Flip/AsyncFlipBegin") and\n');
        fprintf('PsychDataPixx: ERROR: Screen("AsyncFlipEnd/AsyncFlipCheckEnd").\n');
        fprintf('PsychDataPixx: ERROR: This is not allowed! You must call outside that block\n');
        fprintf('PsychDataPixx: ERROR: of code!\n\n');
        error('PsychDataPixx: RequestPsyncedUpdate: Function called too late!');
    end

    % Set flag:
    dpx.needpsync = 1;
    
    return;
end

if strcmpi(cmd, 'IsBusy')
    varargout{1} = dpx.needpsync > 1;
    return;
end

if strcmpi(cmd, 'FlipCount')
    varargout{1} = dpx.swapcount;
    return;
end

if strcmpi(cmd, 'WindowHandle')
    varargout{1} = dpx.window;
    return;
end

if strcmpi(cmd, 'ExecuteAtFlipCount')
    if nargin < 2
        error('PsychDataPixx: ExecuteAtFlipCount: TargetFlipcount argument missing!');
    end
    
    % Assign swapcount of next flip if no target swapcount given, or
    % referring to the past:
    if isempty(varargin{1}) || (varargin{1} <= dpx.swapcount)
        tcount = dpx.swapcount + 1;
    else
        tcount = varargin{1};
    end
    
    if nargin < 3 || isempty(varargin{2})
        error('PsychDataPixx: ExecuteAtFlipCount: CommandString argument missing!');        
    end
    
    % Try to find free slot in schedule, otherwise define a new one, ie.,
    % grow schedule to accomodate additional need:
    freeslot = min(find(dpx.targetswapcounts == -1)); %#ok<MXFND>
    if isempty(freeslot)
        freeslot = length(dpx.targetswapcounts) + 1;
    end
    
    % Add 'tcount' as target swapcount for execution of given command:
    dpx.targetswapcounts(freeslot) = tcount;
    
    % Add given command at corresponding slot of schedule:
    dpx.scheduledcommands{freeslot} = varargin{2};
    
    return;
end

if strcmpi(cmd, 'LoadIdentityClut')
    % Load an identity CLUT into DataPixx at next Screen('Flip'). This is
    % just a little convenience wrapper around 'LoadNormalizedGammaTable':
    if nargin < 2 || isempty(varargin{1})
        error('PsychDataPixx: LoadIdentityClut: Window handle for onscreen window missing!');
    end
    
    linear_lut = repmat(linspace(0, 1, 256)', 1, 3);
    Screen('LoadNormalizedGammaTable', varargin{1}, linear_lut, 2);
    return;
end

if strcmpi(cmd, 'Close')
    % Close of device connection requested.
        
    % Drop refcount:
    dpx.refcount = dpx.refcount - 1;
    
    % No clients left?
    if dpx.refcount <= 0
        % Close device for real:
        doDatapixx('Close');
        
        % Reset whole driver:
        dpx = [];        
    end
    
    return;
end

if strcmpi(cmd, 'Open')
    % Datapixx already open?
    isReady = doDatapixx('IsReady');
    if ~isReady
        % Nope: Open connection to it.
        isReady = doDatapixx('Open');
    end

    % Check if connection is open and ready:
    if isReady
        % Increase our refcount to it:
        dpx.refcount = dpx.refcount + 1;
        
        % Assign default name to device:
        dpx.ID = 'DEFAULT';
        
        % Perform initial syncClocks run:
        syncClocks;
    else
        warning('PsychDataPixx: Device connection not ready after init! TROUBLE AHEAD!!');
    end
    
    return;
end

if strcmpi(cmd, 'SetDummyMode')
    % Set opmode, dummy without device, or the real thing. Apply immediately:
    varargout{1} = dpx.dummymode;
    
    if ~isempty(varargin) && ~isempty(varargin{1})
        dpx.dummymode = varargin{1};
    end
    
    return;
end

if strcmpi(cmd, 'SetVideoMode')
    % Set new video processing mode, apply immediately:
    doDatapixx('SetVideoMode', varargin{1});
    doDatapixx('RegWrRd');
    return;
end

if strcmpi(cmd, 'SetVideoHorizontalSplit')
    % Set new horizontal split mode, apply immediately:
    doDatapixx('SetVideoHorizontalSplit', varargin{1});
    doDatapixx('RegWrRd');
    return;
end

if strcmpi(cmd, 'SetVideoVerticalStereo')
    % Set new vertical sync-doubling mode, apply immediately:
    doDatapixx('SetVideoVerticalStereo', varargin{1});
    doDatapixx('RegWrRd');
    return;
end

% No-Op the new blue-line-stereo and scanning backlight functions, which
% are not yet supported by the Datapixx driver on 64-Bit Octave and
% Matlab on Linux:
if (IsOctave && Is64Bit) || (IsLinux && ~IsOctave)
    % Skip the else branch with the new unsupported commands.
    if strcmpi(cmd, 'EnableVideoStereoBlueline') || strcmpi(cmd, 'DisableVideoStereoBlueline') || ...
       strcmpi(cmd, 'EnableVideoScanningBacklight') || strcmpi(cmd, 'DisableVideoScanningBacklight')
        return;
    end
else
    % New Datapixx driver: Expose new functions...
    
    % Stereo sync line control:
    if strcmpi(cmd, 'EnableVideoStereoBlueline')
        doDatapixx('EnableVideoStereoBlueline');
        doDatapixx('RegWrRd');
        return;
    end
    
    if strcmpi(cmd, 'DisableVideoStereoBlueline')
        doDatapixx('DisableVideoStereoBlueline');
        doDatapixx('RegWrRd');
        return;
    end
    
    % Scanning backlight control:
    if strcmpi(cmd, 'EnableVideoScanningBacklight')
        doDatapixx('EnableVideoScanningBacklight');
        doDatapixx('RegWrRd');
        return;
    end
    
    if strcmpi(cmd, 'DisableVideoScanningBacklight')
        doDatapixx('DisableVideoScanningBacklight');
        doDatapixx('RegWrRd');
        return;
    end
end

if strcmpi(cmd, 'PerformPostWindowOpenSetup')
    % This is called from within PsychImaging()'s PostConfiguration
    % subroutine. It performs all remaining setup work to be done after the
    % onscreen window is open and the imaging pipeline mostly setup:
    
    % Retrieve window handle to onscreen window:
    win = varargin{1};

    if ~isempty(dpx.window)
        error('PsychDataPixx: Tried to open 2nd onscreen window on device while onscreen window already open. There can be only one onscreen window open on the DataPixx display device!');
    end
    
    % Store windowhandle of our window for future use:
    dpx.window = win;
    
    % Load the graphics hardwares gamma table with an identity mapping,
    % so it doesn't interfere with PSYNC code recognition:
    LoadIdentityClut(win);

    % Setup reasonable default PSYNC timeout:
    if isempty(dpx.timeout)
        % Query framerate and calculate a timeout for the device of 5 minutes
        % by default:
        dpx.timeout = Screen('Framerate', win);
        if dpx.timeout == 0
            % OS returns unknown framerate: This means 60 Hz...
            dpx.timeout = 60 * 60 * 5;
        else
            % dpx.timeout fps * 60 secs/minute * 5 minutes:
            dpx.timeout = dpx.timeout * 60 * 5;
        end
    end
    
    % Retrieve window status:
    winfo = Screen('GetWindowInfo', win);
    
    % Retrieve Datapixx status:
    doDatapixx('RegWrRd');
    dpixstatus = doDatapixx('GetVideoStatus');

    % Compute an "all-black" zero pixel line for use in the imagingpipeline
    % callback to fill the first row of the display physical framebuffer
    % with all-zeros in order to facilitate robust operation of the PSYNC
    % mechanism:
    dpx.blackline = uint8(zeros(3, dpixstatus.horizontalResolution));    

    % Tell Datapixx to only accept PSYNC codes in the first scanline (line 0)
    % of the display, and to unconditionally blank that scanline to
    % black, so we avoid distraction visual cues from the PSYNC pixels:
    doDatapixx('SetVideoPixelSyncLine', 0, 1, 1);
    doDatapixx('RegWrRd');
    
    % Attach a callback into the datapixx driver at the end of the
    % left finalizer blit chain to trigger pending preflip operations for
    % DataPixx device:
    Screen('Hookfunction', win, 'AppendMFunction', 'LeftFinalizerBlitChain', 'Final preflip callback into PsychDataPixx driver.', 'PsychDataPixx(0);');
    Screen('HookFunction', win, 'Enable', 'LeftFinalizerBlitChain');

    % Attach a postflip callback for things like timestamp and status
    % collection from the device after flip completion:
    % Screen('Hookfunction', win, 'AppendMFunction', 'ScreenFlipImpliedOperations', 'Postflip callback into PsychDataPixx driver.', 'evalin(''base'', ''PsychDataPixx(2, IMAGINGPIPE_FLIPCOUNT);'');');
    Screen('Hookfunction', win, 'AppendMFunction', 'ScreenFlipImpliedOperations', 'Postflip callback into PsychDataPixx driver.', 'PsychDataPixx(2);');
    Screen('HookFunction', win, 'Enable', 'ScreenFlipImpliedOperations');

    % Attach a window close callback for Device teardown at window close time:
    Screen('Hookfunction', win, 'AppendMFunction', 'CloseOnscreenWindowPostGLShutdown', 'Shutdown window callback into PsychDataPixx driver.', 'PsychDataPixx(''ResetOnWindowClose'');');
    Screen('HookFunction', win, 'Enable', 'CloseOnscreenWindowPostGLShutdown');

    % Stereomode 2 for sync-doubling active?
    if winfo.StereoMode == 2
        % Yes: Need to adapt the y-offset of bottom (right-eye) blitted
        % framebuffer to correct for the extra VBLANK interval which the
        % device inserts in sync doubling mode.

        % Find relevant slot and parameters in imaging pipeline:
        [slot, shaderid, blittercfg, voidptr, glsl, lutid] = Screen('HookFunction', win, 'Query', 'StereoCompositingBlit', 'StereoCompositingShaderCompressedTop');

        % Shader found?
        if slot == -1
            warning('Could not find shader plugin for stereo display inside imaging pipeline for window? Prepare for trouble!!');
        end

        if glsl == 0
            Screen('CloseAll');
            error('Stereo shader is not operational for unknown reason, maybe a Psychtoolbox bug? Sorry...');
        end

        % Build new config string for blitting of left-eye buffer in
        % sync-doubling mode:
        blittercfg = sprintf('Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f:Bilinear', 0, 0, 1.0, 0.5);

        % Insert modified thingy at former position of the old shader:
        Screen('HookFunction', win, 'Remove', 'StereoCompositingBlit', slot);
        posstring = sprintf('InsertAt%iShader', slot);
        Screen('Hookfunction', win, posstring, 'StereoCompositingBlit', shaderid, glsl, blittercfg, lutid);
        
        % Find relevant slot and parameters in imaging pipeline:
        [slot, shaderid, blittercfg, voidptr, glsl, lutid] = Screen('HookFunction', win, 'Query', 'StereoCompositingBlit', 'StereoCompositingShaderCompressedBottom');

        % Shader found?
        if slot == -1
            warning('Could not find shader plugin for stereo display inside imaging pipeline for window? Prepare for trouble!!');
        end

        if glsl == 0
            Screen('CloseAll');
            error('Stereo shader is not operational for unknown reason, maybe a Psychtoolbox bug? Sorry...');
        end

        % Corrective vertical offset for right-eye bottom view: Center it
        % in the middle of the VTOTAL total scanline count, after half of
        % all scanlines. VTOTAL includes the VBL interval of both
        % half-frames:
        voffset = ceil(dpixstatus.verticalTotal / 2) + 0;

        % Build new config string for blitting of right-eye buffer in
        % sync-doubling mode:
        blittercfg = sprintf('Builtin:IdentityBlit:Offset:%i:%i:Scaling:%f:%f:Bilinear', 0, voffset, 1.0, 0.5);

        % Insert modified thingy at former position of the old shader:
        Screen('HookFunction', win, 'Remove', 'StereoCompositingBlit', slot);
        posstring = sprintf('InsertAt%iShader', slot);
        Screen('Hookfunction', win, posstring, 'StereoCompositingBlit', shaderid, glsl, blittercfg, lutid);
    end
    
    return;
end

if strcmpi(cmd, 'CheckGPUSanity')
    
    % Window handle for onscreen window. The imaging pipeline output
    % formatter isn't active at time of this function call:
    win = varargin{1};
    
    % Corrective xoffset for the framebuffer:
    xoffset = varargin{2};

    % Set to 1 for fault injection to test the procedure:
    if 0
        oldlut = Screen('ReadNormalizedGammaTable', win);
        Screen('LoadNormalizedGammaTable', win, oldlut * 0.69);
        %Screen('LoadNormalizedGammaTable', win, rand(size(oldlut)));
    end
    
    % Open our own connection to Datapixx:
    PsychDataPixx('Open');
    
    % Build 256 pixel level RGB test gradient:
    testdata = uint8(zeros(3, 256));
    testdata(:,:) = uint8(repmat(0:255, 3, 1));

    % Clear the onscreen window to black and set its clear color to black.
    Screen('FillRect', win, 0);
    
    % Wait at least two video refresh cycles to make sure the LUT's in the
    % GPU have stabilized:
    doDatapixx('RegWrRdVideoSync');
    doDatapixx('RegWrRdVideoSync');
    
    % At this point the GPU gamma tables should be already loaded with an
    % identity LUT via LoadIdentityClut(). The GL struct and moglcore are
    % initialized so we can use low-level glXXX commands.
    
    % Testing loop: Will display the horizontal test gradient in the top
    % scanline of the display, then read it back from the DataPixx. If
    % everything works correctly then the readback gradient should match
    % the drawn gradient. A mismatch means that either the gamma tables are
    % not loaded with a proper identity LUT, or that some kind of
    % spatio-temporal display dithering is active and interfering.
    %
    % We can try to auto-correct for wrong identity LUT's with our iterative
    % correction loop. It will tweak the gamma tables until a match is
    % achieved. We can't do anything about dithering though beyond
    % reporting it:
    oldlut = Screen('ReadNormalizedGammaTable', win);
    nrlutslots = size(oldlut, 1);
    curlut = oldlut;
    
    psychlasterror('reset');

    % We tweak gamma table elements in steps of 1/1024 units per iteration.
    % As the granularity of Digital DVI-D is only 8 bits or 1/256 units,
    % that means we are resolving 4 times more fine-grained than the device
    % resolution on DVI. This just for safety to avoid "jumping into the
    % wrong bin" due to roundoff errors in the OS/Driver. Could probably do
    % more coarse-grained, but better safe than sorry...
    stepsize = 1/1024;
    nIter = 100;
    failcount = 0;
    successcount = 0;
    regressioncount = 0;
    fluctcount = 0;
    retry = 1;
    
    % Draw test pattern:
    glRasterPos2i(xoffset, 10);
    glDrawPixels(256, 10, GL.RGB, GL.UNSIGNED_BYTE, repmat(testdata, 1, 10));

    % Show it at next retrace:
    Screen('Flip', win);
    fprintf('PsychDataPixx: INFO: GPU gamma table and display encoder test running. Please be patient...\n');
    
    % Enter test and auto-correction loop:
    while retry
        retry2=1;
        while retry2 && (fluctcount < nIter)
            doDatapixx('RegWrRdVideoSync');
            oldrealpixels = doDatapixx('GetVideoLine', 256);
            
            % Repeat measurement until at least 60 consecutive scans give
            % the same value:
            for i=1:60
                doDatapixx('RegWrRdVideoSync');
                realpixels = doDatapixx('GetVideoLine', 256);
                if any(any(int16(oldrealpixels) - int16(realpixels)))
                    % Difference between 1st measurement and i'th
                    % measurement detected. Abort scan and retry:
                    fluctcount = fluctcount + 1;
                    clc;
                    fprintf('PsychDataPixx: INFO: GPU gamma table and display encoder test running. Please be patient...\n');
                    fprintf('Instable DVI-D signal, possibly due to dithering! Will try to fix this by iterative tweaking of gamma tables.\n');
                    fprintf('Iteration Nr. %i in progress. Current pixelvalues received from Datapixx topmost scanline:\n\n', failcount);
                    disp(realpixels);
                    break;
                else
                    if i==60
                        retry2=0;
                    end
                end
            end
        end
        
        if fluctcount >= nIter
            % Totally unstable DVI signal. This could be due to
            % (spatio-)temporal dithering being enabled on the display head
            % or due to serious signal integrity issues.
            fprintf('PsychDataPixx: INFO: No stable signal received during over %i iterations! Interference by dithering?!?\n', nIter);
            fluctcount = 0;     
        end
        
        % 'realpixels' contains a measurement that has been stable for at
        % least 60 video cycles. We use this to judge the quality of the
        % current GPU LUT settings and perform corrections if needed.
        
        % Compute delta matrix to reference values: This should contain
        % all-zeros on properly working GPU's:
        deltavec = int16(testdata) - int16(realpixels);
        
        % Difference?
        if any(any(deltavec))
            % Ohoh!
            if failcount == 0
                fprintf('PsychDataPixx: WARNING: Your GPU has a wrong identity gamma table loaded, possibly due to driver bugs.\n');
                fprintf('PsychDataPixx: WARNING: I will try to automatically compensate for this now. Please standby...\n');
            end
            
            failcount = failcount + 1;
            
            if successcount > 0
                regressioncount = regressioncount + 1;
                fprintf('PsychDataPixx: WARNING: %i. regression detected after tweak %i! Previously successfull setting produces wrong values now\nafter %i successfull iterations!\n', regressioncount, failcount, successcount);
                successcount = 0;
            end
            
            % Convert to double matrix and transpose to match LUT format.
            % Then convert to a corrective increment of stepsize units per
            % 1 pixel unit error:
            deltavec = transpose(double(deltavec)) * stepsize;
            
            % Build correction vector tweakvec: deltavec is a 256 rows by 3
            % columns matrix of delta values to add to curlut for tweaking
            % it to a more correct value. curlut, as read from the GPU,
            % always has 3 columns for (Red, Green, Blue) just like
            % deltavec, but it can have more than 256 slots on high-end
            % GPU's with bigger gamma tables. E.g., the NVidia
            % QuadroFX-3800 when run under Linux with the binary blob
            % driver will have a 2^11 = 2048 slot lut. Therefore we need to
            % "stretch" deltavec vertically so it matches the size of the
            % 'curlut'. By design of the Datapixx and DVI-D data
            % transmission, we can always only compute 256 rows of
            % meaningful tweaking information, therefore we simply repmat
            % replicate successive rows in deltavec to blow it up to a
            % matching tweakvec. This should be good enough to init the 256
            % rows subset of curlut we actually care about with the tweaked
            % values we found:
            tweakvec = zeros(size(curlut));
            chunksize = nrlutslots / 256;
            for tchunk = 0:255
                tweakvec((1 + (tchunk * chunksize)):((tchunk * chunksize) + chunksize), :) = repmat(deltavec(1 + tchunk, :), chunksize, 1);
            end
            
            % Update gamma lut:
            curlut = curlut + tweakvec;
            
            mmcount = sum(sum(deltavec ~= 0));
            fprintf('PsychDataPixx: INFO: Gamma table tweak iteration %i. Had %i mismatching pixels. Retesting with updated LUT...\n', failcount, mmcount);
            
            % Clamp to valid 0-1 range:
            curlut(curlut > 1) = 1;
            curlut(curlut < 0) = 0;
            
            % Upload modified LUT to GPU: We try-catch protect this to
            % avoid error-abort on MS-Windows if invalid gamma tables are
            % encountered:
            try
                Screen('LoadNormalizedGammaTable', win, curlut);
            catch
                psychlasterror('reset');
            end
            
            % Tried too long, too hard?
            if failcount > 255
                % Running with non-zero rasterizer offset?
                if xoffset ~= 0
                    % Hmm. Could be a cascade of bugs. Reset xoffset to
                    % zero and retry whole procedure just to be on the safe
                    % side:
                    fprintf('PsychDataPixx: WARNING: No success so far with rasterizeroffset %i. Retrying with zero-offset...\n', xoffset);
                    failcount = 0;
                    curlut = oldlut;

                    Screen('LoadNormalizedGammaTable', win, oldlut);
                    
                    % Draw test pattern:
                    Screen('Flip', win);
                    xoffset = 0;
                    glRasterPos2i(xoffset, 10);
                    glDrawPixels(256, 10, GL.RGB, GL.UNSIGNED_BYTE, repmat(testdata, 1, 10));

                    % Show it at next retrace:
                    Screen('Flip', win);

                    % Encore une fois...
                    continue;
                end
                
                % Ok, this can't be. Restore gamma table and give up:
                Screen('LoadNormalizedGammaTable', win, oldlut);
                                
                % Fail:
                fprintf('PsychDataPixx: ERROR: Did not reach a stable result after 256 tweak iterations.\n');
                fprintf('PsychDataPixx: ERROR: This is hopeless, something is seriously wrong with your GPU. I give up!\n');
                fprintf('PsychDataPixx: ERROR: One possible reason could be that (spatio-)temporal dithering is enabled for your\n');
                fprintf('PsychDataPixx: ERROR: Datapixx display due to some graphics driver bug, which prevents me from\n');
                fprintf('PsychDataPixx: ERROR: converging to a stable result. Datapixx will not be useable for high-precision\n');
                fprintf('PsychDataPixx: ERROR: color or luminance displays with this broken graphics card/driver. Fix your system!\n');
                
                % Close our own connection to Datapixx:
                PsychDataPixx('Close');
                
                varargout{1} = 1;
                return;
            end
        else
            % No difference. Good. Let's repeat the exercise a few times to
            % rule out bad interference from temporal dithering:
            successcount = successcount + 1;

            fprintf('PsychDataPixx: INFO: Posttest iteration %i passed. Good!\n', successcount)
            
            if successcount > 10
                % 10 successfull consecutive repetitions. That's good
                % enough for us:
                retry = 0;
            end            
        end
    end

    % Ok, we managed to get through:
    Screen('Flip', win);
    
    % How did we do?
    if (successcount > 10) && (failcount == 0)
        % Purrfect! Return full success:
        fprintf('PsychDataPixx: Info: Online test of GPU hardware identity clut and display encoder passed.\n');
        varargout{1} = 0;
    else
        % Got through, but with a non-zero failcount. That probably means
        % that display temporal dithering is not an issue, but that the
        % gamma table needed tweaking.
        fprintf('PsychDataPixx: INFO: Managed to find a useable gamma table setting after %i tweak iterations. :-)\n', failcount);
        if regressioncount > 0
            fprintf('PsychDataPixx: INFO: There were %i regressions throughout the procedure. Could be problems with active display dithering or bad signal quality?\n', regressioncount);
        end
        if fluctcount == 0
            fprintf('\nPsychDataPixx: INFO: There wasn''t any instability during the last testing cycle,\n');
            fprintf('PsychDataPixx: INFO: so harmful display dithering may be absent now and you will be probably fine.\n');
        end
        fprintf('PsychDataPixx: INFO: I will store this LUT on your filesystem for future use.\n');
        if 0
            fprintf('PsychDataPixx: INFO: The difference matrix between the original LUT and the corrected LUT looks like this:\n\n');
            disp(curlut - oldlut);
        end
        fprintf('\n\n');
        
        % Store the current gamma table LUT in a config file for later use
        % by LoadIdentityClut in successive runs:
        SaveIdentityClut(win);
        varargout{1} = 0;        
    end

    % Close our own connection to Datapixx:
    PsychDataPixx('Close');
    
    return;
end

    % Undispatched command? Probably an error:
    error('PsychDataPixx: Unknown command %s! Read the help for valid commands.', cmd);
    return;
end


% Helper function: Call into driver, or perform emulation of real calls in
% dummy mode:
function rc = doDatapixx(varargin)
    global dpx;
    
    if dpx.dummymode
        if (nargin > 0) && strcmpi(varargin{1}, 'GetVideoStatus')
            % Dummy mode: Need to return something for this call: Make up
            % two fake but reasonable values for the returned statusstruct:
            dstat.verticalTotal = 1250;
            dstat.verticalResolution = 1200;
            dstat.horizontalResolution = 1680;
            rc = dstat;
        else
            % Dummy mode: Need to return something for some calls, usually a
            % positive number is the right response, and some want time, so
            % GetSecs always does the right thing ;-)
            rc = GetSecs;
        end
    else
        % Call into MEX-File driver:
        if nargout
            rc = Datapixx(varargin{:});
        else
            rc = [];
            Datapixx(varargin{:});
        end
    end
end

% Clock sync routine: Synchronizes host clock (aka GetSecs time) to box
% internal clock via a sampling and calibration procedure:
function syncresult = syncClocks
    global dpx;
    
    % Default return argument to "invalid":
    syncresult = [nan nan nan];
    
    % Query level of verbosity to use:
    verbosity = dpx.verbosity;
    
    % Perform multiple measurement trials per syncClocks run, take the best
    % one for final result. We use the "best" one because we have a good
    % criterion to find the best one.
    
    % We preallocate the sampling arrays for 250 samples at most. The
    % arrays will grow if 250 samples are not sufficient, at a small
    % performance penalty:
    ntrials = 250;
    
    % Any PSYNC synced write pending?
    if dpx.needpsync
        % Ohoh! Can't syncClocks then!!
        warning('PsychDataPixx: syncClocks: Device is waiting for a Screen(''Flip'') to finish! Can''t sync!');
        return;
    end

    % Dummy mode active?
    if dpx.dummymode
        % Yep, create plausible fake result and return:
        td = GetSecs;
        syncresult = [td td 0.0001];
        dpx.syncSamples(end+1, :) = syncresult;
        return;
    end
    
    % If we reach this point then we are not in emulation mode and can
    % access the real Datapixx MEX-File directly. This is important because
    % call overhead for the emulation may impair the sync accuracy we could
    % get.
    
    % Switch to realtime priority if not already there:
    oldPriority=Priority;
    if oldPriority < MaxPriority('GetSecs')
        Priority(MaxPriority('GetSecs'));
    end
    
    t = zeros(3,ntrials);
    minwin = inf;
    tdeadline = GetSecs + dpx.maxDuration;
    ic = 0;
    
    % Perform measurement trials until success criterion satisfied:
    % Either a sample with a maximum error 'minwin' less than desired
    % threshold, or maximum allowable calibration time reached:
    while (minwin > dpx.optMinwinThreshold) && (GetSecs < tdeadline)
        % Prepare clock query:
        Datapixx('SetMarker');

        % Wait some random fraction of a millisecond. This will desync us
        % from the USB duty cycle and increase the chance of getting a very
        % small time window between scheduling, execution and acknowledge
        % of the send operation:
        WaitSecs(rand / 1000);
        
        % Take pre-Write timestamp in tpre - Sync command not emitted
        % before that time:
        tpre = GetSecs;
        
        % Emit immediate register writecommand:
        Datapixx('RegWr');
        
        % Store completion time in post-write timestamp tpost:
        tpost = GetSecs;
        
        % We know that write command emission has happened at some time
        % after tpre and before tpost. This by design of the USB
        % standard, host controllers and operating system USB stack. This
        % is the only thing we can take for granted wrt. timing, so the
        % "time window" between those two timestamps is our window of
        % uncertainty about the real host time when sync started. However,
        % on a well working system without massive system overload one can
        % be reasonably confident that the real emission of the sync
        % command happened no more than 1 msec before tpost. That is a
        % soft constraint however - useful for computing the final estimate
        % for hosttime, but nothing to be taken 100% for granted.
        
        % Retrieve latched 'SetMarker' timestamp from box, by first reading
        % back the register block, then getting the latched marker value:
        Datapixx('RegWrRd');
        tbox = Datapixx('GetMarker');

        % Compute confidence interval for this sample:
        % For each measurement, the time window tpost - tpre defines a
        % worst case confidence interval for the "real" host system time
        % when the sync command was emitted.
        confidencewindow = tpost - tpre;
        
        % If the confidencewindow is greater than the maximum acceptable
        % window, then we reject this sample, else we accept it:
        if confidencewindow <= dpx.maxMinwinThreshold
            % Within acceptable range. Accept this sample and check if it
            % is the best wrt. to window size so far:
            if confidencewindow < minwin
               % Best confidencewindow so far. Update minwin, as this is one
               % of the abortion critera:
               minwin = confidencewindow;
            end
            
            % Increase sample index to permanently accept this sample for
            % final set of competitors:
            ic = ic + 1;
            
            % Assign values:
            t(1,ic) = tpre;
            t(2,ic) = tpost;
            t(3,ic) = tbox;
        else
            % Inacceptably large error confidencewindow. Reject this sample:
            continue;
        end

        % Next sample pass:
    end

    % Done with sampling: We have up to 'ic' valid samples, unless minwin
    % is still == inf.
    
    % Restore priority to state pre syncClocks:
    if Priority ~= oldPriority
        Priority(oldPriority);
    end
    
    % At least one sample with acceptable precision acquired?
    if (minwin > dpx.maxMinwinThreshold) | (ic < 1) %#ok<OR2>
        % No, not even a single one!
        if verbosity > 1
            fprintf('PsychDataPixx: Warning: On Box "%s", Clock sync failed due to confidence interval of best sample %f secs > allowable maximum %f secs.\n', dpx.ID, minwin, dpx.maxMinwinThreshold);
            fprintf('PsychDataPixx: Warning: Likely your system is massively overloaded or misconfigured!\n');
            fprintf('PsychDataPixx: Warning: See the help for PsychDataPixx(''SyncConstraints'') on how to relax the constraints if you wish to do so.\n');
        end
        
        % That's it:
        return;
    end

    % Ok, we have 'ic' > 0 samples with acceptable precision, according to
    % user specified constraints. Prune result array to valid samples 1 to ic:
    t = t(:, 1:ic);

    % We have three different methods for final sample selection. Each has
    % its favorable cases, but none of them is better by a large margin.
    % All can't guarantee accurate results on highly overloaded systems,
    % but only provide the best under given circumstances.
    %
    % Under low-load conditions, the new style methods may be slightly more
    % accurate under some conditions, but the difference is usually in the
    % sub-millisecond range, so either choice is ok. The options are mostly
    % here for benchmarking and stress testing of the driver...
    %
    % The default method (unless overriden by usercode or
    % userconfiguration) is method 1 -- New style with postwrite timestamps.
    % This because this method empirically seems to give best results on
    % low or medium load on our test set of machines and operating systems.
    switch dpx.syncmode
        case 2
            % Old style method 2 - Middle of minwin confidence window:
            
            % Select sample with smallest confidence interval [tpre; tpost]:
            [mintdiff, idx] = min(t(2,:) - t(1,:));

            % Host time corresponds to midpoint of the confidence interval,
            % assuming a uniform distribution of likelyhood of the true write
            % time in the interval [tpre ; tpost]:
            hosttime = (t(1,idx) + t(2,idx)) / 2;
            
        case 1
            % New style method 1 - Postwrite timestamps:

            % Choose the most accurate sample from the set of candidates. This is
            % the sample with the smallest difference between the postwrite
            % timestamp and the associated box timestamp, ie., with the smallest
            % offset between postwrite host clock time and box clock time at
            % receive of sync command. The reasoning behind this goes like this:
            %
            % 1) The time offset between host clock and box clock is a constant -
            % at least within a significant multi-second time interval between
            % successive syncClocks calls (due to only small clock drift), but
            % certainly within a syncClocks run of a few hundred milliseconds
            % (error due to clock drift in this interval is negligible).
            %
            % 2) Computed box clock time t(3,:) is "close to perfect", as this
            % timestamp is taken by box microprocessor and firmware with a very
            % small and basically constant delay after sync token receive, ie.,
            % write completion. (Maximum theoretical error is smaller than 0.1 msecs).
            %
            % 3) The correct and optimal clock offset between host and box would be
            % tdiff = tsend - t(3,:) iff tsend would be host time at true write
            % completion.
            %
            % 4) The measured host time at write completion t(2,:) is always later
            % (and therefore numerically greater) than the true host time tsend at
            % write completion due to an unknown, random, greater than zero delay
            % tdelta, i.e., t(2,:) = tsend + tdelta, tdelta > 0. tdelta is the sum
            % of:
            %
            % a) The unknown delay of up to 1 msec between USB write-URB completion
            % by the USB host controller (which would be the real completion time
            % tsend) and detection of completion due to USB IOC (Interrupt-On-
            % Completion) due to invocation of the host controllers hardware
            % interrupt handler and host controller schedule scan and URB
            % retirement inside the interrupt handler.
            %
            % b) Random (and theoretically unbounded) scheduling delay / execution
            % delay between status update of the serial port data structures by the
            % interrupt handler and detection of write completion + timestamping by
            % the IOPort driver in polling mode, or scheduling delay between
            % wakeup-operation caused by the interrupt handler and start of
            % execution of the timestamping in the IOPort driver in blocking mode.
            %
            % The syncClocks error is therefore directly proportional to the size
            % of tdelta. Therefore:
            %
            % tdiff(:) = t(2,:) - t(3,:) by definition of clock offset host vs. box.
            % t(2,:) = tsend(:) + tdelta(:) by unknown scheduling/execution noise tdelta.
            %
            % It follows that by above definitions:
            %
            % tdiff(:) = tsend(:) - t(3,:) + tdelta(:);
            %
            % --> As we defined tsend(:) to be the unknown, but perfect and
            % noise-free, true send timestamp, and t(3,:) to be the perfect receive
            % timestamp by the box, it follows that by selecting the sample 'idx'
            % with the minimal tdiff(idx) from the set tdiff(:), we will select the
            % sample with the unknown, but minimal tdelta(idx). As tdelta accounts
            % for all the remaining calibration error, minimizing tdelta will
            % maximize the accuracy of the clock sync.
            %
            % ==> Select sample with minimum t(2,:) - t(3,:) as final best result:
            [mintdiff, idx] = min(t(2,:) - t(3,:));

            % mintdiff is our best estimate of clock offset host vs. box, and
            % t(:,idx) is the associated best sample. Unfortunately there isn't any
            % way to compute the exact residual calibration error tdelta(idx). The
            % only thing we know is that the error is bounded by the length of the
            % associated 'minwin' confidence interval of this sample, so we will
            % return 'minwin' as an upper bound on the calibration error. As
            % 'minwin' was used as a threshold in the sample loop for outlier
            % rejection, we can be certain that our estimate carries no greater
            % error than 'dpx.maxMinwinThreshold'.

            % Extract all relevant values for the final sample:

            % Host time corresponds to tpost write timestamp, which should be as
            % close as possible to real host send timestamp:
            hosttime = t(2,idx);

        case 0
            % New style method 0 - Prewrite timestamps:

            % Choose the most accurate sample from the set of candidates.
            % This is the sample with the biggest difference between the
            % prewrite timestamp and the associated box timestamp, ie.,
            % with the biggest offset between prewrite host clock time and
            % box clock time at receive of sync command. The reasoning
            % behind this goes like this:
            %
            % 1) The time offset between host clock and box clock is a constant -
            % at least within a significant multi-second time interval between
            % successive syncClocks calls (due to only small clock drift), but
            % certainly within a syncClocks run of a few hundred milliseconds
            % (error due to clock drift in this interval is negligible).
            %
            % 2) Computed box clock time t(3,:) is "close to perfect", as this
            % timestamp is taken by box microprocessor and firmware with a very
            % small and basically constant delay after sync token receive, ie.,
            % write completion. (Maximum theoretical error is smaller than 0.1 msecs).
            %
            % 3) The correct and optimal clock offset between host and box would be
            % tdiff = tsend - t(3,:) where tsend would be unknown host time
            % at true write completion.
            %
            % 4) The measured host time before write submission t(1,:) is
            % always earlier
            % (and therefore numerically smaller) than the true host time tsend at
            % write completion due to an unknown, random, greater than zero delay
            % tdelta, i.e., t(1,:) = tsend - tdelta, tdelta > 0. tdelta is the sum
            % of:
            %
            % a) The unknown delay of up to 1 msec between USB write-URB completion
            % by the USB host controller (which would be the real completion time
            % tsend) and detection of completion due to USB IOC (Interrupt-On-
            % Completion) due to invocation of the host controllers hardware
            % interrupt handler and host controller schedule scan and URB
            % retirement inside the interrupt handler.
            %
            % b) Random (and theoretically unbounded) scheduling delay / execution
            % delay between status update of the serial port data structures by the
            % interrupt handler and detection of write completion + timestamping by
            % the IOPort driver in polling mode, or scheduling delay between
            % wakeup-operation caused by the interrupt handler and start of
            % execution of the timestamping in the IOPort driver in blocking mode.
            %
            % The syncClocks error is therefore directly proportional to the size
            % of tdelta. Therefore:
            %
            % tdiff(:) = t(1,:) - t(3,:) by definition of clock offset host vs. box.
            % t(1,:) = tsend(:) - tdelta(:) by unknown execution noise tdelta.
            %
            % It follows that by above definitions:
            %
            % tdiff(:) = tsend(:) - t(3,:) - tdelta(:);
            %
            % --> As we defined tsend(:) to be the unknown, but perfect and
            % noise-free, true send timestamp, and t(3,:) to be the perfect receive
            % timestamp by the box, it follows that by selecting the sample 'idx'
            % with the maximal tdiff(idx) from the set tdiff(:), we will select the
            % sample with the unknown, but minimal tdelta(idx). As tdelta accounts
            % for all the remaining calibration error, minimizing tdelta will
            % maximize the accuracy of the clock sync.
            %
            % ==> Select sample with maximum t(1,:) - t(3,:) as final best result:
            [maxtdiff, idx] = max(t(1,:) - t(3,:));

            % maxtdiff is our best estimate of clock offset host vs. box, and
            % t(:,idx) is the associated best sample. Unfortunately there isn't any
            % way to compute the exact residual calibration error tdelta(idx). The
            % only thing we know is that the error is bounded by the length of the
            % associated 'minwin' confidence interval of this sample, so we will
            % return 'minwin' as an upper bound on the calibration error. As
            % 'minwin' was used as a threshold in the sample loop for outlier
            % rejection, we can be certain that our estimate carries no greater
            % error than 'dpx.maxMinwinThreshold'.

            % Extract all relevant values for the final sample:

            % Host time corresponds to tpre write timestamp, which should be as
            % close as possible to real host send timestamp:
            hosttime = t(1,idx);
        otherwise
            error('PsychDataPixx: syncClocks: Unknown timestamping method provided. This is a driverbug!!');
    end
    
    % Box timers time taken "as is":
    boxtime  = t(3,idx);
    
    % Recalculate upper bound on worst case error 'minwin' from this best
    % samples tpost - tpre:
    minwin = t(2,idx) - t(1,idx);

    if verbosity > 3
        fprintf('PsychDataPixx: ClockSync(%i) on box "%s": Got %i valid samples, maxconfidence interval = %f msecs, winner interval %f msecs.\n', dpx.syncmode, dpx.ID, ic, 1000 * dpx.maxMinwinThreshold, 1000 * minwin);
        if dpx.syncmode == 1
            fprintf('PsychDataPixx: Confidence windows in interval [%f - %f] msecs. Range of clock offset variation: %f msecs.\n', 1000 * min(t(2,:)-t(1,:)), 1000 * max(t(2,:)-t(1,:)), 1000 * range(t(2,:) - t(3,:)));
        else
            fprintf('PsychDataPixx: Confidence windows in interval [%f - %f] msecs. Range of clock offset variation: %f msecs.\n', 1000 * min(t(2,:)-t(1,:)), 1000 * max(t(2,:)-t(1,:)), 1000 * range(t(1,:) - t(3,:)));
        end
    end
    
    % Assign (host,box,confidence) sample to sync struct:
    syncresult = [hosttime, boxtime, minwin];
    
    % Also store the sample in the syncSamples history:
    dpx.syncSamples(end+1, :) = syncresult;
    
    % SyncClocks run finished:
end

% Map box timestamps to host clock (aka GetSecs) timestamps, based on clock
% sync results from syncClocks, and clock drift correction from clock ratio
% calibration:
function timing = box2GetSecsTime(timing)
    global dpx;

    thostbase = dpx.syncSamples(end, 1);
    tboxbase  = dpx.syncSamples(end, 2);

    % Compute delta in box seconds from time when calibration was
    % done, wrt. box timer, then multiply with clockRatio
    % conversion factor to translate that delta into a delta wrt.
    % host clock, so timing will be the delta in host clock seconds
    % wrt. calibration time:
    timing = (timing - tboxbase) * dpx.clkRatio;

    % Add basetime on host - this should give us correctly mapped
    % time in the GetSecs timebase of the computer:
    timing = timing + thostbase;
end

% Map box timestamps to host clock (aka GetSecs) timestamps, based on clock
% sync results from all syncClocks samples from a session:
function [timing, sd, clockratio] = box2GetSecsTimePostHoc(timing)
    global dpx;
    
    % Check if the latest syncClocks sample is older than 30 seconds. If
    % so, then we acquire a new final sample. We also resample if the last
    % sample is of too low accuracy, or if there are less than 2 samples in
    % total, as the fitting procedure needs at least 2 samples to work:
    while (size(dpx.syncSamples, 1) < 2) || ...
          ((GetSecs - dpx.syncSamples(end, 1)) > 30) || ...
          (dpx.syncSamples(end, 3) > dpx.maxMinwinThreshold)

      % Perform a syncClocks to get a fresh sample to finalize the sampleset:
      syncClocks;
    end

    % Extract samples for line fit:
    tbox  = dpx.syncSamples(:, 2);
    thost = dpx.syncSamples(:, 1);
    
    % MK: Change of strategy: We always use polyfit() instead of sometimes
    % regress() or robustfit() for the following reasons:
    %
    % 1. Strategy switching may introduce artifacts / non-predictable
    % behaviour into the results.
    %
    % 2. A simple polyfit() is better understood than the more advanced
    % robustfit() with its tons of tunable parameters and strategies, for
    % which i don't know which are appropriate for our data.
    %
    % 3. robustfit() and regress() require the Matlab statistics toolbox
    % and therefore won't work on a plain standard Matlab installation,
    % whereas polyfit() is part of default Matlab and Octave.
    %
    % The old code is left here for documentation purpose:
    %
    %     % More than 2 samples available?
    %     if 0 && length(tbox) > 2
    %         % Octave and older versions of Matlab don't have 'robustfit',
    %         % so we fall back to 'regress' if this function is lacking:
    %         if exist('robustfit') %#ok<EXIST>
    %             [coef st]=robustfit(tbox,thost);  % fit a line
    %             sd=st.robust_s; % stddev. in seconds.
    %         else
    %             coef =regress(thost, [ones(size(thost,1), 1), tbox ]);  % fit a line
    %             sd=0; % stddev. undefined with regress().
    %         end
    %         clockratio = coef(2);
    %         clockbias  = coef(1);
    %     else
    
    % Use polyfit, i.e., a least squares fit of a polynomial of either
    % degree 1 (purely linear), or degree 2 (quadratic). We choose a purely
    % linear fit if only 2 samples are available - the best we can do. If
    % more than 2 samples are available, we use a quadratic fit to account
    % for small non-linearities in the clock hardware over long intervals:
    [coef, st, mu] = polyfit(tbox, thost, min(length(tbox)-1, 2));

    clockratio = coef(1);
    clockbias  = coef(2); %#ok<NASGU>

    sd.norm = st.normr;

    %     end
    
    % Ok, got mapping equation getsecst = timing * clockratio + clockbias;
    % Apply it to our input timestamps:
    %
    %    timing = timing * clockratio + clockbias;
    %
    %
    % Actually, we use polyval() to evaluate the equation. For linear case,
    % this is the same as above equation, but it conveniently also handles
    % the quadratic case and provides some margins on mapping error:
    if length(tbox) > 3
        [timing, sd.delta] = polyval(coef, timing, st, mu);
    else
        % sd.delta undefined for less than 4 elements in timing:
        sd.delta = repmat(Inf, size(timing));
        timing = polyval(coef, timing, st, mu);
    end

    if dpx.verbosity > 3
        fprintf('PsychDataPixx: In post-hoc box->host mapping: Coefficients of mapping polynom (in decreasing order) are: ');
        disp(coef);
        fprintf('\n\n');
    end
    
    % Ready.
end
