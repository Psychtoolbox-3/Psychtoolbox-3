function PsychPortAudioDataPixxTimingTest(waitTime, exactstart, deviceid, latbias, triggerLevel)
% PsychPortAudioDataPixxTimingTest([waitTime = 1][, exactstart=1][, deviceid=-1][, latbias=0][, triggerLevel=0.01])
%
% Test script for sound onset timing reliability and sound onset latency of
% the PsychPortAudio sound driver.
%
% This script configures the driver for low latency and high timing
% precision, then executes ten trials where it tries to emit a beep sound,
% starting in exact sync to a black-white transition on the display.
%
% You'll need measurement equipment to use this: A DataPixx device from
% VPixx technologies, connected to your computer via the USB connection
% cable. Also a connection between the line-out jack of your soundcard and
% the line-in jack of the DataPixx to transmit the sound data. The DataPixx
% will receive the audio output of PsychPortAudio/Your Soundcard, timestamp
% it and send the computed timing data to your computer via USB.
%
% Some parameters may need tweaking, especially the 'triggerLevel'
% parameter.
%
% Optional parameters:
%
% 'waitTime'   = Time to wait (in seconds) before playing sound. Defaults
%                to 1 second if omitted.
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
% 'triggerLevel' = Sound signal amplitude for DataPixx to detect sound
%                  onset. Defaults to 0.01 = 1% of max amplitude if
%                  exactstart == 0, otherwise it is auto-detected by
%                  calibration. This will likely need tweaking on your
%                  setup.
%

nTrials = 10;

% Initialize driver, request low-latency preinit:
InitializePsychSound(1);

% Force GetSecs and WaitSecs into memory to avoid latency later on:
GetSecs; WaitSecs(0.1);

if nargin < 1
    waitTime = [];
end

if isempty(waitTime)
    waitTime = 1;
end

% If 'exactstart' wasn't provided, assume user wants to test exact sync of
% audio and video onset, instead of testing total onset latency:
if nargin < 2
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
if nargin < 3
    deviceid = [];
end

if isempty(deviceid)
    deviceid = -1;
end

% Needs to determined via measurement once for each piece of audio
% hardware:
if nargin < 4
    latbias = [];
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

% Request latency mode 2:
reqlatencyclass = 2;

% Requested output frequency, may need adaptation on some audio-hw:
freq = 48000;       % Must set this. 96khz, 48khz, 44.1khz.
buffersize = 0;     % Pointless to set this. Auto-selected to be optimal.
suggestedLatencySecs = [];

if IsARM
    % ARM processor, probably the RaspberryPi SoC. This can not quite handle the
    % low latency settings of a Intel PC, so be more lenient:
    suggestedLatencySecs = 0.025;
    if isempty(latbias)
        latbias = 0.000593;
        fprintf('Choosing a latbias setting of 0.000593 secs or 0.593 msecs, assuming this is a RaspberryPi ARM SoC.\n');
    end
    fprintf('Choosing a high suggestedLatencySecs setting of 25 msecs to account for lower performing ARM SoC.\n');
end

if isempty(latbias)
    % Unknown system: Assume zero bias. User can override with measured
    % values:
    fprintf('No "latbias" provided. Assuming zero bias. You''ll need to determine this via measurement for best results...\n');
    latbias = 0;
end

if nargin < 5
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

% Perform one warmup trial, to get the sound hardware fully up and running,
% performing whatever lazy initialization only happens at real first use.
% This "useless" warmup will allow for lower latency for start of playback
% during actual use of the audio driver in the real trials:
PsychPortAudio('Start', pahandle, 1, 0, 1);
PsychPortAudio('Stop', pahandle, 1);

% Ok, now the audio hardware is fully initialized and our driver is on
% hot-standby, ready to start playback of any sound with minimal latency.

% Switch to realtime scheduling at maximum allowable Priority:
Priority(MaxPriority(0));

% Initialize audio capture subsystem of Datapixx:
% 96 KhZ sampling rate, Mono capture: Average across channels (0), Audio
% input is line-in (2), Gain is 1.0 (1):
DatapixxAudioKey('Open', 96000, 0, 2, 1);

% Check settings by printing them:
dpixstatus = Datapixx('GetMicrophoneStatus') %#ok<NOPRT,NASGU>

% Auto-Selection of triggerLevel for Datapixx timestamping requested?
if exactstart && isempty(triggerLevel)
    % Use auto-trigger mode. Tell the function how long the silence
    % interval at start of each trial is expected to be. This will be
    % used for calibration: We set it to 75% of the duration of the pause
    % between start of Datapixx recording and scheduled sound onset time:
    DatapixxAudioKey('AutoTriggerLevel', waitTime * 0.75);
    fprintf('Setting lead time of silence in Datapixx auto-trigger mode to %f msecs.\n', waitTime * 0.75 * 1000);    
else
    % Triggerlevel for DataPixx sound onset detection:
    if isempty(triggerLevel)
        % Default to 1%:
        triggerLevel = 0.01;
    end

    fprintf('Using a trigger level for DataPixx of %f. This may need tweaking by you...\n', triggerLevel);
    DatapixxAudioKey('TriggerLevel', triggerLevel);
end

% Wait for keypress.
fprintf('\n\nPress any key to start measurement.\n\n');
KbStrokeWait;

% nTrials measurement trials:
for i=1:nTrials
    % Start audio capture on DataPixx now. Return true 'tStartBox'
    % timestamp of start in box clock time:
    tStartBox = DatapixxAudioKey('CaptureNow');

    if exactstart
        % Schedule start of audio at exactly 'waitTime' seconds ahead:
        PsychPortAudio('Start', pahandle, 1, GetSecs + waitTime, 0);
        desired = waitTime;
    else
        % No test of scheduling, but of absolute latency: Start audio
        % playback immediately:
        PsychPortAudio('Start', pahandle, 1, 0, 0);
        desired = 0;
    end

    if 0
        % Spin-Wait until hw reports the first sample is played...
        offset = 0; %#ok<*UNRCH>
        while offset == 0
            status = PsychPortAudio('GetStatus', pahandle);
            offset = status.PositionSecs;
            plat = status.PredictedLatency;
            fprintf('Predicted Latency: %6.6f msecs.\n', plat*1000);
            if offset>0
                break;
            end
            WaitSecs('YieldSecs', 0.001);
        end
    end

    % Retrieve true delay from DataPixx measurement and stop recording on the device:
    [audiodata, measuredAudioDelta] = DatapixxAudioKey('GetResponse', waitTime + 1, [], 1); %#ok<*ASGLU>

    % Compute expected delay based on audio onset time as predicted/measured by
    % PsychPortAudio:
    status = PsychPortAudio('GetStatus', pahandle);
    tPortAudio(i) = status.StartTime; %#ok<AGROW>
    tDataPixx(i)  = tStartBox + measuredAudioDelta; %#ok<AGROW>

    fprintf('Buffersize %i, xruns = %i, playpos = %6.6f secs, measured audio onset delay = %6.1f msecs vs. desired %6.1f\n', status.BufferSize, status.XRuns, status.PositionSecs, 1000 * measuredAudioDelta, 1000 * desired);

    if 0
        figure;
        plot(audiodata);
    end

    % Stop playback:
    PsychPortAudio('Stop', pahandle, 1);
end

% Remap Datapixx clock timestamps to Psychtoolbox GetSecs() timestamps:
tDataPixx = PsychDataPixx('BoxsecsToGetsecs', tDataPixx);

% Done: Back to normal scheduling:
Priority(0);

% Close Datapixx audio subsystem:
DatapixxAudioKey('Close');

% Close PsychPortAudio:
PsychPortAudio('Close');

fprintf('\n\n');
for i =1:nTrials
    audioDelta(i) = 1000 * (tDataPixx(i) - tPortAudio(i)); %#ok<AGROW>
    fprintf('%i. PsychPortAudio measured onset timestamp error is %6.6f msecs.\n', i, audioDelta(i));
end

% Discard 1st trial:
audioDelta = audioDelta(2:end);
fprintf('\nAvg timestamp error %6.6f msecs, Stddev %6.6f msecs, Range %6.6f msecs.\n\n', mean(audioDelta), std(audioDelta), range(audioDelta));

return;
