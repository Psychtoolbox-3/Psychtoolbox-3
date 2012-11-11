function PsychPortAudioTimingTest(exactstart, deviceid, latbias, waitframes, useDPixx, triggerLevel)
% PsychPortAudioTimingTest([exactstart=1][, deviceid=-1][, latbias=0][, waitframes=1][, useDPixx=0][, triggerLevel=0.01])
%
% Test script for sound onset timing reliability and sound onset
% latency of the PsychPortAudio sound driver.
%
% This script configures the driver for low latency and high timing
% precision, then executes ten trials where it tries to emit a beep sound,
% starting in exact sync to a black-white transition on the display.
%
% You'll need measurement equipment to use this: A photo-diode attached to
% the top-left corner of your CRT display, a microphone attached to your
% speakers, some oszillograph to record and measure the signals from the
% diode and microphone.
%
% Some parameters may need tweaking. Make sure you got a setup as described
% in 'help InitializePsychSound' for best results.
%
% Optional parameters:
%
% 'exactstart' = 0 -- Start immediately, measure absolute latency.
%              = 1 -- Test accuracy of scheduled sound onset. (Default)
%
% 'deviceid'   = -1 -- Auto-select optimal device (Default).
%             >=0   -- Select specified output device. See
%                      PsychPortAudio('GetDevices') for a list of devices.
%
% 'latbias'    = Hardware inherent latency bias. To be determined by
%                measurement - allows to PA to correct for it if provided.
%                Unit is seconds. Defaults to zero.
%
% 'waitframes' = Time to wait (in video refresh intervals) before emitting beep + flash.
%                Defaults to a value of 1, may need to be set higher for
%                high latency sound hardware (where absolute latency > 1
%                video refresh duration).
%
% 'useDPixx'   = 1 -- Use DataPixx device to automatically measure the true
%                     audio onset time wrt. to visual stimulus onset.
%                0 -- Don't use DataPixx. This is the default.
%
% 'triggerLevel' = Sound signal amplitude for DataPixx to detect sound
%                  onset. Defaults to 0.01 = 1% of max amplitude if
%                  exactstart == 0, otherwise it is auto-detected by
%                  calibration. This will likely need tweaking on your
%                  setup. If the measured audio onset delta by DataPixx is
%                  much lower (or almost zero) than the expected delta
%                  reported by PsychPortAudio, then the triggerLevel may be
%                  too low and you should try if slightly higher thresholds
%                  help to discriminate signal from noise. Too high values
%                  may cause a hang of the script. In practice, levels
%                  between 0.01 and 0.1 should yield good results. Setting
%                  the 'useDPixx' flag to 2 also plots the waveforms
%                  captured by DataPixx, which may help in selection of the
%                  optimal triggerLevel.
%

% Initialize driver, request low-latency preinit:
InitializePsychSound(1);

if ~IsLinux
  PsychPortAudio('Verbosity', 10);
end

% Force GetSecs and WaitSecs into memory to avoid latency later on:
GetSecs;
WaitSecs(0.1);

% If 'exactstart' wasn't provided, assume user wants to test exact sync of
% audio and video onset, instead of testing total onset latency:
if nargin < 1 
   exactstart = [];
end

if isempty(exactstart)
    exactstart = 1;
end

if exactstart
    fprintf('Will test accuracy of scheduled sound onset, i.e. how well the driver manages to\n');
    fprintf('emit sound at exactly the specified "when" deadline. Sound should start in exact\n');
    fprintf('sync with display black-white transition (or at least very close - < 1 msec off).\n');
    fprintf('The remaining bias can be corrected by providing the bias as "latbias" parameter\n');
    fprintf('to this script. Variance of sound onset between trials should be very low, much\n');
    fprintf('smaller than 1 millisecond on a well working system.\n\n');
else
    fprintf('Well test total latency for immediate start of sound playback, i.e., the "when"\n');
    fprintf('parameter is set to zero. The difference between display black-white transition\n');
    fprintf('and start of emitted sound will be the total system latency.\n\n');
end

% Default to auto-selected default output device if none specified:
if nargin < 2 
   deviceid = [];
end

if isempty(deviceid)
    deviceid = -1;
end

if deviceid == -1
    fprintf('Will use auto-selected default output device. This is the system default output\n');
    fprintf('device in "normal" (=reliable but high latency) mode. In low-latency mode its the\n');
    fprintf('device with the lowest inherent latency on your system (as determined by some internal\n');
    fprintf('heuristic). If you are not satisfied with the results you may query the available devices\n');
    fprintf('yourself via a call to devs = PsychPortAudio(''GetDevices''); and provide the device index\n');
    fprintf('of a suitable device\n\n');
else
    fprintf('Selected the following output device (deviceid=%i) according to your spec:\n', deviceid);
    devs = PsychPortAudio('GetDevices');
    for idx = 1:length(devs)
        if devs(idx).DeviceIndex == deviceid
            break;
        end
    end
    disp(devs(idx));
end

% Request latency mode 2, which used to be the best one in our measurement:
% classes 3 and 4 didn't yield any improvements, sometimes they even caused
% problems.
reqlatencyclass = 2;

% Requested output frequency, may need adaptation on some audio-hw:
freq = 44100;       % Must set this. 96khz, 48khz, 44.1khz.
buffersize = 0;     % Pointless to set this. Auto-selected to be optimal.
suggestedLatencySecs = [];

if IsWin
    % Hack to accomodate bad Windows systems or sound cards. By default,
    % the more aggressive default setting of something like 5 msecs can
    % cause sound artifacts on cheaper / less pro sound cards:
    suggestedLatencySecs = 0.015 %#ok<NOPRT>
    fprintf('Choosing a high suggestedLatencySecs setting of 15 msecs to account for shoddy Windows operating system.\n');
    fprintf('For low-latency applications, you may want to tweak this to lower values if your system works better than average timing-wise.\n');
end

% Needs to determined via measurement once for each piece of audio
% hardware:
if nargin < 3 
   latbias = [];
end

if isempty(latbias)
    % Unknown system: Assume zero bias. User can override with measured
    % values:
    fprintf('No "latbias" provided. Assuming zero bias. You''ll need to determine this via measurement for best results...\n');
    latbias = 0;
end

if nargin < 4
    waitframes = [];
end

if nargin < 5
    useDPixx = [];
end

if isempty(useDPixx)
    useDPixx = 0;
end

if nargin < 6
    % Default triggerLevel is "auto-trigger":
    triggerLevel = [];
end

% Open audio device for low-latency output:
pahandle = PsychPortAudio('Open', deviceid, [], reqlatencyclass, freq, 2, buffersize, suggestedLatencySecs);

% Tell driver about hardwares inherent latency, determined via calibration
% once:
prelat = PsychPortAudio('LatencyBias', pahandle, latbias) %#ok<NOPRT,NASGU>
postlat = PsychPortAudio('LatencyBias', pahandle) %#ok<NOPRT,NASGU>

% Generate some beep sound 1000 Hz, 0.1 secs, 50% amplitude:
mynoise(1,:) = 0.5 * MakeBeep(1000, 0.1, freq);
mynoise(2,:) = mynoise(1,:);

% Fill buffer with data:
PsychPortAudio('FillBuffer', pahandle, mynoise);

% Setup display:
screenid = max(Screen('Screens'));

% Shall we use the DataPixx for measurement?
if useDPixx
    % Yes!
    
    % Initialize audio capture subsystem of Datapixx:
    % 96 KhZ sampling rate, Mono capture: Average across channels (0), Audio
    % input is line in (2), Gain is 1.0 (1):
    DatapixxAudioKey('Open', 96000, 0, 2, 1);
    
    % Check settings by printing them:
    dpixstatus = Datapixx('GetMicrophoneStatus') %#ok<NOPRT,NASGU>

    if ~(exactstart && isempty(triggerLevel))
        if isempty(triggerLevel)
            % Choose a default of 1% of max. signal amplitude:
            triggerLevel = 0.01;
        end
        fprintf('Using a trigger level for DataPixx of %f. This may need tweaking by you...\n', triggerLevel);
        DatapixxAudioKey('TriggerLevel', triggerLevel);
    end
    
    % DataPixx: Setup Screen imagingpipeline to support measurement via the PSYNC
    % video synchronization mode of DataPixx and Screen():
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseDataPixx');
    win = PsychImaging('OpenWindow', screenid, 0);
    LoadIdentityClut(win);
else
    % Default: No need for imaging pipeline:
    win = Screen('OpenWindow', screenid, 0);
end

ifi = Screen('GetFlipInterval', win);

% Set waitframes to a good default, if none is provided by user:
if isempty(waitframes)
    % We try to choose a waitframes that maximizes the chance of hitting
    % the onset deadline. We are conservative in our estimate, because a
    % few video refresh cycles hardly matter for this test, but increase
    % our chance of success without need for manual tuning by user:
    if isempty(suggestedLatencySecs)
        % Let's assume 12 msecs on Linux and OSX as a achievable latency by
        % default, then double it:
        waitframes = ceil((2 * 0.012) / ifi) + 1;        
    else
        % Whatever was provided, then double it:
        waitframes = ceil((2 * suggestedLatencySecs) / ifi) + 1;
    end
end

fprintf('\n\nWaiting %i video refresh cycles before white-flash.\n', waitframes);

% Auto-Selection of triggerLevel for Datapixx timestamping requested?
if useDPixx && exactstart && isempty(triggerLevel)
    % Use auto-trigger mode. Tell the function how long the silence
    % interval at start of each trial is expected to be. This will be
    % used for calibration: We set it to 75% of the duration of the pause
    % between start of Datapixx recording and scheduled sound onset time:
    DatapixxAudioKey('AutoTriggerLevel', ifi * waitframes * 0.75);
    fprintf('Setting lead time of silence in Datapixx auto-trigger mode to %f msecs.\n', ifi * waitframes * 0.75 * 1000);
end

% Perform one warmup trial, to get the sound hardware fully up and running,
% performing whatever lazy initialization only happens at real first use.
% This "useless" warmup will allow for lower latency for start of playback
% during actual use of the audio driver in the real trials:
PsychPortAudio('Start', pahandle, 1, 0, 1);
PsychPortAudio('Stop', pahandle, 1);

% Ok, now the audio hardware is fully initialized and our driver is on
% hot-standby, ready to start playback of any sound with minimal latency.

% Wait for keypress.
KbStrokeWait;

% Realtime scheduling: Can be used if otherwise timing is not good enough.
% Priority(MaxPriority(win));

% Ten measurement trials:
for i=1:10
    if useDPixx
        % Schedule start of audio capture on DataPixx at next Screen('Flip'):
        DatapixxAudioKey('CaptureAtFlip');
    end
    
    % This flip clears the display to black and returns timestamp of black onset:
    % It also triggers start of audio recording by the DataPixx, if it is
    % used, so the DataPixx gets some lead-time before actual audio onset.
    [vbl1 visonset1]= Screen('Flip', win);

    % Prepare black white transition:
    Screen('FillRect', win, 255);
    Screen('DrawingFinished', win);

    if exactstart
        % Schedule start of audio at exactly the predicted visual stimulus
        % onset caused by the next flip command.
        PsychPortAudio('Start', pahandle, 1, visonset1 + waitframes * ifi, 0);
    end

    % Ok, the next flip will do a black-white transition...
    [vbl visual_onset t1] = Screen('Flip', win, vbl1 + (waitframes - 0.5) * ifi);

    if ~exactstart
        % No test of scheduling, but of absolute latency: Start audio
        % playback immediately:
        PsychPortAudio('Start', pahandle, 1, 0, 0);
    end

    t2 = GetSecs;

    % Spin-Wait until hw reports the first sample is played...
    offset = 0;
    while offset == 0
        status = PsychPortAudio('GetStatus', pahandle);
        offset = status.PositionSecs;
        t3=GetSecs;
        plat = status.PredictedLatency;
        fprintf('Predicted Latency: %6.6f msecs.\n', plat*1000);
        if offset>0
            break;
        end
        WaitSecs('YieldSecs', 0.001);
    end
    audio_onset = status.StartTime;

    %fprintf('Expected visual onset at %6.6f secs.\n', visual_onset);
    %fprintf('Sound started between %6.6f and  %6.6f\n', t1, t2);
    %fprintf('Expected latency sound - visual = %6.6f\n', t2 - visual_onset);
    %fprintf('First sound buffer played at %6.6f\n', t3);
    fprintf('Flip delay = %6.6f secs.  Flipend vs. VBL %6.6f\n', vbl - vbl1, t1-vbl);
    fprintf('Delay start vs. played: %6.6f secs, offset %f\n', t3 - t2, offset);

    fprintf('Buffersize %i, xruns = %i, playpos = %6.6f secs.\n', status.BufferSize, status.XRuns, status.PositionSecs);
    fprintf('Screen    expects visual onset at %6.6f secs.\n', visual_onset);
    fprintf('PortAudio expects audio onset  at %6.6f secs.\n', audio_onset);
    fprintf('Expected audio-visual delay    is %6.6f msecs.\n', (audio_onset - visual_onset)*1000.0);

    if useDPixx
        % 'visonset1' is the GetSecs() time of start of capture on
        % DataPixx. 'audio_onset' is reported GetSecs() audio onset time
        % according to PsychPortAudio.
        %
        % 'expectedAudioDelta' is therefore the expected delay for the
        % measured audio onset by DataPixx:
        expectedAudioDelta = audio_onset - visonset1;
        
        % Retrieve true delay from DataPixx measurement and stop recording
        % on the device:
        [audiodata, measuredAudioDelta] = DatapixxAudioKey('GetResponse', [], [], 1);
        fprintf('DPixx: Expected audio onset delta is %6.6f secs.\n', expectedAudioDelta);
        fprintf('DPixx: Measured audio onset delta is %6.6f secs.\n', measuredAudioDelta);
        fprintf('DPixx: PsychPortAudio measured onset error is therefore %6.6f msecs.\n', 1000 * (measuredAudioDelta - expectedAudioDelta));
        if useDPixx > 1
            figure;
            plot(audiodata);
        end
    end
    
    % Stop playback:
    PsychPortAudio('Stop', pahandle, 1);

    % Wait a bit...
    WaitSecs(0.3);

    Screen('FillRect', win, 0);
    telapsed = Screen('Flip', win) - visual_onset; %#ok<NASGU>
    WaitSecs(0.6);
end

% Done, close driver and display:
Priority(0);

if useDPixx
    % Close Datapixx audio subsystem:
    DatapixxAudioKey('Close');
end

PsychPortAudio('Close');
Screen('CloseAll');

% Done. Bye.
return;
