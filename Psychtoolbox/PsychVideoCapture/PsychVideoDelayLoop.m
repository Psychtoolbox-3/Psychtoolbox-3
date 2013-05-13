function rc = PsychVideoDelayLoop(cmd, varargin)
% PsychVideoDelayLoop(subcommand, arg1, arg2, ...)
%
% This implements a realtime video feedback loop with adjustable
% delay, e.g., for action-perception studies.
%
% Arguments:
% subcommand - Is a string containing the subcommand to call.
% arg1, ... argn - Are the arguments that the specific subcommand
% requires.
%
% Subcommands, their meaning and arguments:
%
% PsychVideoDelayLoop('Verbosity', level);
% -- Set level of verbosity: 0 == Shut up. 1 == Errors and warnings.
% 2 == Information as well. The default is 1.
%
% handle = PsychVideoDelayLoop('Open', windowPtr [,deviceId] [,ROI] [,inColor])
% -- Open a video capture device 'deviceId' for use on a specific onscreen
% window 'windowPtr', setup region of interest 'ROI', select if capture
% should happen 'inColor' = 1 or gray-scale (inColor = 0).
% This returns a 'handle' for the device, so other scripts can do
% device specific setup.
%
% PsychVideoDelayLoop('Close')
% -- Shutdown, close and release all video capture devices. This
% invalidates any handle obtained from 'Open'.
%
% fps = PsychVideoDelayLoop('TuneVideoRefresh', capturerate)
% -- Measure the exact capture rate of the video device when
% requesting a capturerate of 'capturerate', check if the displays
% refresh rate is roughly compatible with the capture device setting
% and then try to fine-tune the display refresh rate in order to
% minimize phase-shifts between capture device work-cycle and display
% device work-cycle. Return the final fine-tuned and measured framerate
% 'fps' of the display device.
% Note: Fine-tuning is only possible in a very narrow range (+/- 1Hz)
% around the display refresh rate set in your display settings. This
% feature is currently only supported on GNU/Linux.
%
% PsychVideoDelayLoop('SetAbortKeys', keyarray)
% -- Define a sequence of keycodes for valid abort
% keys. If any of the keys given in the sequence is pressed, the
% video loop will exit. You can map keys to keycodes via KbName, e.g.,
% key1 = KbName('Escape'); key2 = KbName('Space'); keyarray = [key1 key2];
% --> Create a keyarray that would abort on Escape- or Space- keypress.
%
% PsychVideoDelayLoop('SetAbortTimeout', timeout)
% -- Define a maximum duration of the feedback loop in seconds. After
% that amount of time has elapsed, the loop will exit.
%
% PsychVideoDelayLoop('SetHeadstart', timemargin)
% -- Define an estimate of how long the system will take to process a
% new video frame plus the expected drift during one trial. The loop
% will take this extra amount of time into account when prestarting the
% camera to make sure that the excess latency caused by the video loop itself
% is as short and stable as possible. E.g., a value of 0.004 secs for
% processing overhead + maybe 0.004 secs for drift == 0.008 secs could
% be reasonable for an otherwise well synchronized system.
%
% PsychVideoDelayLoop('SetPresentation', fullfov, mirrored, upsidedown);
% -- Change mode of presentation: fullfov=0 Show centered image, fullfov=1
% Zoom image to fill full area of onscreen window. mirrored=0 Normal
% presentation, mirrored=1 Mirror left-right. upsidedown=0 Upright,
% upsidedown=1 Upside-Down.
%
% PsychVideoDelayLoop('SetLogging', mode, maxseconds)
% -- Disable (mode=0) or enable (mode=1) logging of timestamps.
% 'maxseconds' is the number of seconds for which the arrays should
% be pre-allocated. After running the delay loop you can query the
% logged timestamps via the 'GetLog' subfunction.
%
% log = PsychVideoDelayLoop('GetLog')
% -- Return the timing logs of last loop run. This is a 3 rows by
% n columns matrix, where each column corresponds to the timing
% samples of one frame: log(1,i) = Absolute system time in seconds,
% when frame i was captured. log(2,i) = Delta (seconds) between
% capture and visual onset of image i on screen. log(3,i) contains
% an estimate of the full delay between capture onset of frame i and
% visual onset. It is the value log(2,i) + estimated latency between
% start of camera sensor exposure and transmit completion for the frame.
% The accuracy of this estimate should be pretty good for cameras known
% to Psychtoolbox and it is a *guess* for the lower bound on the real
% latency on unknown cameras.
%
% PsychVideoDelayLoop('RecordFrames', framestep)
% -- Record every 'framestep'th frame in system RAM as an OpenGL
% texture. The video loop doesn't discard every texture after drawing
% it, but enqueues it an a system RAM buffer. The vector of texture
% handles can be retrieved via 'GetRecordedFrames' after the loop
% has finished. Default is zero == Recording disabled. A value of
% one records every frame, a value of two every second, ...
%
% Calling this function will also reset the vector of recorded frames
% to empty, but it will not delete the textures! That is your task.
%
% texids = PsychVideoDelayLoop('GetRecordedFrames')
% Return vector of texture handles for all recorded frames.
% texids(1,i) contains the texture handle for the i'th recorded frame.
% texids(2,i) contains the onset timestamp for the i'th recorded frame.
%
% PsychVideoDelayLoop('RunLoop', delayFrames[, onlinecontrol]);
% -- Run the video feedback loop, using the parameters specified above.
% The video loop will start and run until one of the abort keys is pressed,
% or the timeout is reached. It will log timestamps as requested. Each
% captured video frame is output again after 'delayFrames' capture cycle
% durations, e.g., capturerate = 30 fps -> cycle = 1/30 sec = 33.33 ms -->
% delay is at least delayFrames * 33.33 ms. Images are drawn and shown
% in sync with vertical retrace after that amount of time. The real onset
% time obviously depends on the monitor refresh interval and phase between
% camera and monitor.
%
% If 'onlinecontrol' is set to 1, then a few control keys are enabled to
% allow for interactive change of settings like brightness, gain and
% exposure time: 'b' increases brightness, 'd' decreases brightness.
% Up/DownArrow keys increase/decrease gain. Right/LeftArrow keys
% increase/decrease exposure time (shutter time).
%

% History:
% 8.08.06  Written (MK)
% 15.02.10 Small improvements and fixes (MK)

% Window handle of target window:
persistent win;
persistent ifi;
persistent fps;
persistent grabber;
persistent capturerate;
persistent cfi;
persistent ROI;
persistent camlatency;
persistent headstart;
persistent recordframes;
persistent captexs;
persistent isOpen;
persistent abortkeys;
persistent aborttimeout;
persistent fullfov;
persistent mirrored;
persistent upsidedown;
persistent maxlogduration;
persistent logmode;
persistent timestamps;
persistent exposureinc;
persistent exposuredec;
persistent gaininc;
persistent gaindec;
persistent brightinc;
persistent brightdec;
persistent verbose;

if isempty(isOpen)
    isOpen = 0;
    KbName('UnifyKeyNames');
    exposureinc = KbName('RightArrow');
    exposuredec = KbName('LeftArrow');
    gaininc = KbName('UpArrow');
    gaindec = KbName('DownArrow');
    brightinc = KbName('b');
    brightdec = KbName('d');
end

if isempty(verbose)
    verbose = 1;
end;

if isempty(aborttimeout)
    aborttimeout = Inf;
end;

if isempty(fullfov)
    fullfov = 1;
end;

if isempty(mirrored)
    mirrored = 0;
end;

if isempty(upsidedown)
    upsidedown = 0;
end;

if isempty(logmode)
    logmode = 0;
end;

if isempty(maxlogduration)
    maxlogduration = 1000;
end;

if isempty(camlatency)
    camlatency = 0;
end;

if isempty(headstart)
    headstart = 0.005 + 0.0042;
end;

if isempty(recordframes)
    recordframes = -1;
    captexs=[];
end;

if nargin < 1
    error('Subcommand missing: You need to spec at least one subcommand!');
end;

if strcmp(cmd, 'Verbosity')
    if nargin < 2
        error('You must provide a level of verbosity when calling subcommand ''Verbosity'' ');
    end

    verbose = varargin{1};
    rc = 0;
    return;
end

if strcmp(cmd, 'Open')
    % Subcommand 'Open': Open video capture device.
    if isOpen
        error('Called ''Open'' although video capture device already opened!');
    end

    if nargin < 2 || isempty(varargin(1))
        error('In Open: You must specify a valid onscreen window handle.');
    end;
    win = varargin{1};

    if nargin < 3 || isempty(varargin(2))
        deviceId = [];
    else
        deviceId = varargin{2};
    end;

    if nargin < 4 || isempty(varargin(3))
        ROI = [0 0 640 480];
    else
        ROI = varargin{3};
    end;

    if nargin < 5 || isempty(varargin(4))
        pixelsize = 1;
    else
        if varargin{4}>0
            pixelsize = 3;
        else
            pixelsize = 1;
        end;
    end;

    % Query refresh interval of display:
    ifi = Screen('GetFlipInterval', win);

    % Open video capture device with given parameters. We want the default
    % number of buffers, but we disallow use of the slow fallback-path.
    grabber = Screen('OpenVideoCapture', win, deviceId, ROI, pixelsize, [], 0);
    rc = grabber;

    isOpen = 1;

    % triggercount = Screen('SetVideoCaptureParameter', grabber, 'WaitTriggerCount')

    % Done with opening the device. Return devicehandle in rc:
    return;
end

if strcmp(cmd, 'Close')
    % Close capture device:
    if ~isOpen
        error('You need to call ''Open'' on a display and video capture device first!');
    end

    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);

    isOpen = 0;

    return;
end

if strcmp(cmd, 'SetAbortKeys')
    if nargin < 2
        error('You need to specify the array with keycodes for abort keys!');
    end

    abortkeys = varargin{1};

    return;
end

if strcmp(cmd, 'SetAbortTimeout')
    if nargin < 2
        error('You need to specify the timeout (in seconds) for Runloop!');
    end

    aborttimeout = varargin{1};

    if aborttimeout <= 0
        error('You need to specify a positive (greater than zero) timeout for Runloop!');
    end

    return;
end

if strcmp(cmd, 'SetHeadstart')
    if nargin < 2
        error('You need to specify the time margin (in seconds) in SetHeadstart!');
    end

    headstart = varargin{1};

    if headstart < 0
        error('You need to specify a positive (greater or equal zero) time for SetHeadstart!');
    end

    return;
end

if strcmp(cmd, 'RecordFrames')
    if nargin < 2
        error('You need to specify the framestep argument in RecordFrames!');
    end

    recordframes = round(varargin{1});
    captexs = [];

    if recordframes <= 0
        % A value of minus one signals recording disabled. This value is
        % crucial for some optimizations in the video loop to reduce
        % loop execution overhead.
        recordframes = -1;
    end

    return;
end

if strcmp(cmd, 'GetRecordedFrames')
    rc = captexs;
    return;
end

if strcmp(cmd, 'SetPresentation')
    if nargin < 4
        error('You need to specify all three arguments for SetPresentation!');
    end

    fullfov = varargin{1};
    mirrored = varargin{2};
    upsidedown = varargin{3};

    return;
end

if strcmp(cmd, 'SetLogging')
    if nargin < 2
        error('You need to specify at least the logging mode for SetLogging!');
    end

    logmode = varargin{1};

    if nargin >= 3
        maxlogduration = varargin{2};
    end;

    return;
end

if strcmp(cmd, 'GetLog')
    rc = timestamps;
    return;
end

if strcmp(cmd, 'TuneVideoRefresh')
    % Tuning of video refresh rate of display to capture device requested.
    if ~isOpen
        error('You need to call ''Open'' on a display and video capture device first!');
    end

    if nargin < 2
        error('You need to specify the requested capture frame rate for the video capture device!');
    end

    capturerate = varargin{1};

    % Start capture on our capture device:
    Screen('StartVideoCapture', grabber, capturerate, 1);

    % Throw away the first 2 frames:
    Screen('GetCapturedImage', win, grabber, 2);
    Screen('GetCapturedImage', win, grabber, 2);
    [dummy oldcts] = Screen('GetCapturedImage', win, grabber, 2);

    % Measurement loop: Run 100 frames:
    cfi = 0;
    for i = 1:100
        % We just retrieve the capture timestamp, nothing else:
        [dummy cts] = Screen('GetCapturedImage', win, grabber, 2);
        delta = cts - oldcts;
        oldcts = cts;
        cfi = cfi + delta;
    end

    % Estimate camera latency from current capture settings.
    camlatency = PsychCamSettings('EstimateLatency', grabber);

    Screen('StopVideoCapture', grabber);

    % cfi is the measured interval between two captured images:
    cfi = cfi / 100;
    capturerate = 1/cfi;

    % We need a display refresh rate that is roughly a multiple of 'capturerate'
    fps = Screen('Framerate', win, 1);
    ratemultiplier = round(fps / capturerate);

    % targetfps would be the optimal display refresh rate for our camera:
    targetfps = ratemultiplier * capturerate;

    % Try to set display to 'targetfps':
    realfps = Screen('Framerate', win, 2, targetfps);

    if realfps<=0
        % Failed to set requested rate. Try to get as close as possible:
        if targetfps > fps
            % Direction -1 --> Need to get faster.
            direction = -1;
        else
            % Direction +1 --> Need to slow down.
            direction = +1;
        end

        realfps = 1;
        while realfps>0
            realfps = Screen('Framerate', win, 2, direction);
        end
    end

    % Ok, we are as close as possible: Query result.
    realfps = Screen('Framerate', win, 1);

    % Good enough?
    if abs(realfps - targetfps)>1
        % Nope :(
        if verbose>0
            fprintf('Failed to achieve a good match between capture framerate and display framerate.\n');
            fprintf('Try to manually set your display to a refresh rate which is a multiple of %f !\n', capturerate);
        end
        % Restore old refresh rate:
        Screen('Framerate', win, 2, fps);
    else
        % Yes.
        if verbose>1
            fprintf('Display fps changed from %f Hz to %f Hz to match good targetfps of %f Hz. Residual mismatch %f Hz.\n', fps, realfps, targetfps, realfps - targetfps);
        end
    end

    % Recalibrate display: Take 100 valid samples:
    ifi = Screen('GetFlipInterval', win, 100);
    fps = 1/ifi;

    capturerate = varargin{1};

    % Return new fps value:
    rc = fps;
    return;
end

if strcmp(cmd, 'RunLoop')
    % Run the video delay loop.
    if ~isOpen
        error('You need to call ''Open'' on a display and video capture device first!');
    end

    if nargin < 2
        error('delayFrames parameter missing for RunLoop!');
    end

    % Get fifodelay in units of captureslots:
    fifodelay = varargin{1};

    if fifodelay < 0
        error('You need to specify a positive delayFrames parameter for RunLoop!');
    end

    if nargin >= 3
        onlinecontrol = varargin{2};
    else
        onlinecontrol = 0;
    end

    % Reset keycode to zero:
    keycode = zeros(1, 256);

    % Query real monitor refresh interval (in seconds) as computed by
    % Psychtoolbox internal calibration:
    ifi = Screen('GetFlipInterval', win);

    % Recompute estimated camera latency, just to be safe...
    camlatency = PsychCamSettings('EstimateLatency', grabber, capturerate);

    % Translate latency in frames into latency in milliseconds:
    latencymillisecs = fifodelay * cfi * 1000.0;
    if verbose > 1
        fprintf('Requested minimum delay in milliseconds: %f\n', latencymillisecs);
        fprintf('Estimated camera latency with current settings is %f ms.\n', camlatency * 1000);
    end;

    % Allocate timestamp buffers for 'maxlogduration' seconds:
    if logmode > 0
        timestamps= zeros(3, min(maxlogduration, aborttimeout+1) * capturerate);
    end

    % Allocate texture id recording vector, if recording enabled.
    captexscount = 0;
    if recordframes > 0
        captexs = zeros(2, ceil((aborttimeout+1) * capturerate / recordframes));
    end

    % Allocate ringbuffer of texture handles and capture timestamps
    % for implementation of the delay loop:
    videofifo = zeros(2, fifodelay+1);

    % Preinit capture counter: It runs 'fifodelay' frames ahead of read counter:
    capturecount = fifodelay;

    % Setup our recycled texture handle: Initially it is zero which means
    % "No old unused texture available for recycling"...
    recycledtex = 0;

    % Select a stepwidth for exposuretime of 1 for unknown cams, 0.1 ms for
    % known cams.
    expdelta=1;
    if PsychCamSettings('IsKnownCamera', grabber)
        expdelta = 0.1;
    end

    % Monoscopic view (windowed or fullscreen):
    if fullfov>0
        dstrect = Screen('Rect', win);
    else
        dstrect = CenterRect(ROI, Screen('Rect', win));
    end

    if mirrored>0 || upsidedown>0
        % Setup transformation for our image.
        Screen('glPushMatrix', win);
        [hw hh] = RectCenter(dstrect);

        Screen('glTranslate', win, hw, hh, 0);
        if mirrored>0
            sx = -1;
        else
            sx = 1;
        end;

        if upsidedown>0
            sy = -1;
        else
            sy = 1;
        end;

        Screen('glScale', win, sx, sy, 1);
        Screen('glTranslate', win, -hw, -hh, 0);

    end;

    % Fetch and throw away all stale frames that are pending in the grabber FIFO:
    texid = 1;
    while texid > 0
        texid = Screen('GetCapturedImage', win, grabber, 0);
        if texid>0
            Screen('Close', texid);
        end;
    end;

    % Sync us to the retrace and get a retrace timestamp:
    synctime = Screen('Flip', win) + (ceil(1.0 / ifi) * ifi);

    % Compute start time to be 1 second from now, minus some slack that
    % is necessary to compensate for camera latency and processing delay,
    % so our first frame will arrive as close to retrace as possible, but
    % with some security margin to account for system drift and scheduling
    % jitter.
    synctime = synctime - headstart - camlatency;

    % Start video capture with 'capturerate' frames per second, if possible. Use
    % low-latency capture by dropping frames, if necessary. Try to start at system time 'synctime':
    [capturerate starttime] = Screen('StartVideoCapture', grabber, capturerate, 1, synctime);
    % Startdelta is difference between real start time and requested start time:
    startdelta = starttime - synctime;

    % Record start time of feedback loop and sync us to VBL:
    % We do a double-flip, just to clear out both framebuffers.
    Screen('Flip', win);
    tstart=Screen('Flip', win);
    oldcts = tstart;
    tonset = tstart;

    % Video capture and feedback loop. Runs until keypress or 'aborttimeout' secs have passed:
    while (tonset - tstart) < aborttimeout
        % Calling KbCheck is expensive (often more than 1 millisecond). We avoid it,
        % if onlinecontrol is disabled and no abort keys are set.
        if ~isempty(abortkeys) || onlinecontrol
            % Check for keypress:
            [down secs keycode] = KbCheck;
            
            if down
                % Key pressed. Check which one and process it:
                
                % Any of the abort-keys pressed?
                if any(intersect(find(keycode), abortkeys))
                    % Abort key pressed: Exit the loop.
                    break;
                end;

                % None of the abort-keys pressed. Online control enabled?
                if onlinecontrol>0
                    % Yes. Check if a control key is pressed and handle it:
                    if keycode(exposureinc)
                        value = PsychCamSettings('ExposureTime', grabber) + expdelta;
                        PsychCamSettings('ExposureTime', grabber, value);
                        camlatency = PsychCamSettings('EstimateLatency', grabber);

                    end

                    if keycode(exposuredec)
                        value = PsychCamSettings('ExposureTime', grabber) - expdelta;
                        PsychCamSettings('ExposureTime', grabber, value);
                        camlatency = PsychCamSettings('EstimateLatency', grabber);
                    end

                    if keycode(gaininc)
                        value = PsychCamSettings('Gain', grabber) + 1;
                        PsychCamSettings('Gain', grabber, value);
                    end

                    if keycode(gaindec)
                        value = PsychCamSettings('Gain', grabber) - 1;
                        PsychCamSettings('Gain', grabber, value);
                    end

                    if keycode(brightinc)
                        value = PsychCamSettings('Brightness', grabber) + 1;
                        PsychCamSettings('Brightness', grabber, value);
                    end

                    if keycode(brightdec)
                        value = PsychCamSettings('Brightness', grabber) - 1;
                        PsychCamSettings('Brightness', grabber, value);
                    end

                    if verbose > 1
                        fprintf('Estimated camera latency with current settings is %f ms.\n', camlatency * 1000);
                    end
                end;
            end;
        end; % Of keyboard checking and processing...

        %        mytelapsed1 = GetSecs - tonset

        % Retrieve most recently captured image from video source, block if none is
        % available yet. If recycledtex is a valid handle to an old, no longer needed
        % texture, the capture engine will recycle it for higher efficiency:
        [tex cts nrdropped]=Screen('GetCapturedImage', win, grabber, 1, recycledtex);

        % Frame dropped during this capture cycle?
        if nrdropped > 0 && verbose > 1
            fprintf('Frame dropped: %i, %f\n', capturecount, (cts-oldcts)*1000);
        end;

        if (cts < oldcts && oldcts~=tstart && verbose >0)
            fprintf('BUG! TIMESTAMP REVERSION AT FRAME %i, DELTA = %f !!!\n', capturecount, (cts-oldcts)*1000);
        end;

        oldcts=cts;

        % New image captured and returned as texture?
        if (tex>0)
            % Yes. Put it into our fifo ringbuffer, together with the requested presentation
            % deadline for that image, which is the capture timestamp + requested delay:
            capturecount = capturecount + 1;
            writeptr = privateMod(capturecount, size(videofifo, 2)) + 1;
            videofifo(1, writeptr) = tex;
            videofifo(2, writeptr) = cts + (latencymillisecs / 1000.0);
            if logmode>0
                % Store capture timestamp in seconds of system time.
                timestamps(1, capturecount)=cts;
            end

            % Done with capture of this frame...

            % recycledtex has been used up. Null it out:
            recycledtex = 0;

            % Now read out the frame some fifo delayslots behind and show it:
            readcount = capturecount - fifodelay;
            readptr = privateMod(readcount, size(videofifo, 2)) + 1;

            % Get texture handle for image to show:
            tex = videofifo(1, readptr);

            % Null-out this used up texture in video fifo.
            videofifo(1, readptr)=0;

            % Nothing to show yet?
            if tex == 0
                % Skip remaining loop:
                continue;
            end;

            % Draw the image.
            Screen('DrawTexture', win, tex, [], dstrect);

            %          mytelapsed2 = Screen('DrawingFinished', win, 2, 1)

            % Perform image onset in sync with retrace and get onset timestamp.
            tonset = Screen('Flip', win, 0, 2);
            if logmode > 0
                % Compute and log the delay between capture and display.
                timestamps(2, readcount) = tonset - timestamps(1, readcount);
                timestamps(3, readcount) = timestamps(2, readcount) + camlatency;
            end;

            if (recordframes>0) && (privateMod(readcount, recordframes)==0)
                captexscount = captexscount + 1;
                captexs(1, captexscount) = tex;
                captexs(2, captexscount) = tonset;
            else
                % We do not need texture 'tex' anymore. Put it into 'recycledtex',
                % so the framecapture engine can reuse it for faster processing.
                recycledtex = tex;
            end;
        end;

        % We have processed the fifo content. Repeat the loop to see
        % if new frames are ready for display.
    end;

    % Final flip. Clears the backbuffer to background color:
    Screen('Flip', win);

    % Done with video feedback loop. Shutdown video capture:
    rc.telapsed = GetSecs - tstart;

    if mirrored>0 || upsidedown>0
        % Undo image transforms:
        Screen('glPopMatrix', win);
    end;

    % Stop capture, do the stats:
    rc.droppedincapturedevice = Screen('StopVideoCapture', grabber);

    rc.avgfps = readcount / rc.telapsed;

    % Truncate vector of texture indices:
    if captexscount>=1
        captexs = captexs(:, 1:captexscount);
    end;

    % Set invalid entries at beginning of timestamps to zero:
    if fifodelay>0
        timestamps(2:3, 1:fifodelay)=0;
    end;

    % Truncate timestamps array to its real size:
    if size(timestamps,2) > readcount
        timestamps = timestamps(:, 1:readcount);
    end;

    rc.keycode = keycode;
    rc.totaldisplayed = readcount;
    rc.startdelta = startdelta;

    % Well done!
    return;
end

help PsychVideoDelayLoop;

fprintf('\n\nUnknown subcommand: %s\n', cmd);
error('Unknown subcommand specified! Please read help text above.');

end

function rem = privateMod(x,y)
% Our private modulo implementation. Only handles positive scalar
% values correctly, but should be much faster than Matlabs/Octaves
% general mod(x,y) function.
rem = x - (y * floor(x / y));
end
