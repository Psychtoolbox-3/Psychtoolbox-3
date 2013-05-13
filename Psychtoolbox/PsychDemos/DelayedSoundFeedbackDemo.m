function DelayedSoundFeedbackDemo(reqlatency, duplex, freq, minLatency)
% DelayedSoundFeedbackDemo([reqlatency=150 ms][, duplex=0][, freq = 48000][, minLatency= 10 ms])
%
% THIS IS EARLY ALPHA CODE! IT MAY OR MAY NOT WORK RELIABLY ON YOUR SETUP!
% TEST IT WITH MEASUREMENT EQUIPMENT IF YOU DEPEND ON ACCURATE FEEDBACK
% TIMING!!!
%
% Demonstrates usage of the new Psychtoolbox sound driver PsychPortAudio()
% for audio feedback with a controlled delay.
%
% Sound is captured from the default recording device and then - with a
% selectable delay - played back via the default output device.
%
% Parameters and their meaning:
%
% 'reqlatency' Wanted feedback latency between sound input and output in
% milliseconds. A value of zero will ask for the lowest possible latency on
% the given setup. Defaults to 150 msecs. Please notice that the minimum
% achievable latency will be constrained by the capabilities of your
% operating system, sound card driver, computer hardware and sound
% hardware. Only very high quality systems will be able to go below 5 msecs
% latency, good systems will be able to go below 20 msecs, but less capable
% setups may only allow for a latency much larger than 20 msecs. In order
% to achieve low latency reliably without timing glitches or audible
% artifacts, you may need to tune both the parameters for this demo and
% your system setup carefully. The optimal parameter set varies from setup
% to setup.
%
% 'duplex' = Select between full-duplex and half-duplex mode:
%
% Depending on your sound hardware you'll have to either leave 'duplex' at
% its default of zero (2 times half-duplex mode, aka simplex mode) or set
% it to 1 (full-duplex mode). On a given system, only one of these will work
% reliably (or at all): ASIO audio hardware -- typically on MS-Windows --
% will usually need full-duplex mode and won't work at all in simplex mode.
% On Macintosh OS/X it depends on the sound hardware. IntelMacs are happy
% with half-duplex mode, some PowerMacs may need full-duplex mode. However,
% except for 'reqlatency' == 0 minimal latency mode, simplex mode provides
% much higher accuracy and reliability on OS/X at least with the built-in
% soundchips on Intel based Macintosh computers. On Linux, performance
% varies depending on the card at use.
%
% 'freq' = Sampling frequency (Hz). Defaults to 48000 Hz as this rate is
% most commonly supported on most sound hardware. The maximum achievable
% value depends on your specific soundcard. IntelMac's built in soundchips
% allow for a maximum of 96000 Hz, high-end soundcards may allow for 192000
% Hz under some circumstances. Increasing the frequency reduces minimum
% latency but increases system load and the probability of glitches.
%
% 'minLatency' is a tuning parameter for the driver and a hard-constraint
% on the mininum achievable latency for feedback. It is ignored on OS/X,
% but can be tused for tuning latency vs. reliability on Linux and on
% MS-Windows. High-end cards may allow for much lower than the default 10
% msecs, low-end cards may malfunction at lower settings. Non-ASIO
% soundcards on MS-Windows will likely fail already at much higher settings
% and be therefore completely unsuitabe for low latency feedback.
%
% Specific tips for different setups:
%
% On OS/X with builtin soundchip on IntelMacs, choose duplex = 0 for
% feedback with controlled low-latency, and a freq'ency of 96000 Hz. For
% lowest latency mode, you may try reqlatency = 0 and duplex = 1.
%
% On MS-Windows you *must* use a soundcard with ASIO support for any
% reasonable results! Here you should always set duplex = 1 for full-duplex
% operation, anything else will fail. Use reqlatency = 0 for feedback with
% minimal latency, positive values for feedback with controlled latency.
% Play around with the 'minLatency' parameter, set it as low as possible -
% to the lowest value that doesn't cause any error messages by our driver
% or audible artifacts like crackling noises or static. Try to set
% 'freq'uency as high as possible. Check the manual of your soundcard for
% the highest value that can be used for capture + playback. E.g., the
% Soundblaster Audigy ZS 2 seems to be limited to max. 48000 Hz in this
% mode.
%
% On Linux, no general statements can be made, only that some soundcards
% allow for extremely low latencies of < 2 msecs if properly tuned. Search
% the Internet for tips.
%

%
% If you need low-latency, make sure to read "help InitializePsychSound"
% carefully or contact the forum.
%

% History:
% 07/20/2007 Written (MK)
% 07/19/2009 Derived and largely rewritten from BasicSoundFeedbackDemo (MK)
% 08/02/2009 Add support for 'DirectInputMonitoring' for reqLatency=0 (MK)
% 04/03/2011 Disable dynamic adaptation of captureQuantum. It causes
%            artifacts at some feedback delay settings. (MK)

% Level of debug output:
verbose = 1;

% Close plots from previous invocation, if any:
close all;

% Store diagnostic timestamps in a huge array:
tc = 0;
tstats = zeros(4, 1000000);

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

% Preload GetSecs, need it later in time-critical part:
GetSecs;

fprintf('\n\nTHIS IS EARLY ALPHA CODE! IT MAY OR MAY NOT WORK RELIABLY ON YOUR SETUP!\nTEST IT WITH MEASUREMENT EQUIPMENT IF YOU DEPEND ON ACCURATE FEEDBACK\nTIMING!!!\n\n');

% Only check ESCape key in KbCheck to save some hazzle and computation time:
KbName('UnifyKeyNames');
RestrictKeysForKbCheck(KbName('ESCAPE'));

% Latency provided?
if nargin < 1
    reqlatency = [];
end

if isempty(reqlatency)
    reqlatency = 150;
end

if nargin < 2
    duplex = [];
end

if isempty(duplex)
    duplex =0;
end

if nargin < 3
    freq = [];
end

if isempty(freq)
    % Default to a sample rate of 48kHz. Should be supported by most hardware:
    freq = 48000;
end

if nargin < 4
    minLatency = [];
end

if isempty(minLatency)
    % Default to a safety margin of 10 msecs if no other provided: It won't
    % be possible to achieve latencies less than minLatency, regardless
    % what you do:
    minLatency = 10;
end

% Map requested minLatency from msecs to seconds:
minLatency = minLatency / 1000;

if ~IsWin
    % Set minLatency to [] aka auto-select on any system but Windows:
    % The PsychPortAudio driver knows best what to select on the other OS,
    % whereas there's no good rule for Windows + ASIO, so user must do
    % trial and error until she finds a value that is low enough for low
    % latency, but high enough for stability and glitch-free audio:
    minLatency = [];
end

% Map requested latency from msecs to seconds:
lat =  reqlatency / 1000;

% Wait for release of all keys on keyboard:
KbReleaseWait;

% Perform low-level initialization of the sound driver:
InitializePsychSound(1);

% Select the most aggressive latency mode latmode = 4. This will try to get
% low-latency and exact timing at all costs, even if it disrupts or crashes
% other running sound applications and consumes lots of ressources. It will
% fail with an error if it can't achieve the most high-perf settings:
latmode = 4;

% Provide some debug output:
PsychPortAudio('Verbosity', 10);

if (reqlatency == 0) && duplex
    % Special case: Full-duplex mode with minimum latency. We bypass Matlab
    % by activating PsychPortAudios full-duplex monitoring mode. The driver
    % itself will feed back all captured sound to the outputs with lowest
    % possible latency. However we don't have any control over latency or
    % sound and this only works on full-duplex hardware...
    %
    % More specifically: It should work well on Windows + ASIO hardware,
    % but not very well on at least OS/X:
    pa = PsychPortAudio('Open', [], 4+2+1, latmode, freq, 2, [], minLatency);

    % Now that the device is open, try to enable the "Zero latency direct input
    % monitoring" feature of some subset of ASIO V2.0 compliant cards. If
    % supported, this feature will cause the card to route sound directly
    % from its input connectors to its output connectors without any
    % extended processing on the card itself, and without any processing by
    % the host computer or our driver. A typical implementation on a
    % higher-end card would configure the analog mixing and amplifier
    % circuits on the card to directly route analog from input to output,
    % without even AD/DA conversion. This would provide true zero latency
    % feedback.
    %
    % This call will ask to enable (=1) routing of all inputs (inputChannel
    % = -1) to output 0 and following (outputChannel = 0), applying zero
    % decibel gain (gain = 0.0) ie., no attenuation or amplification of the
    % signal. On stereo outputs, signals shall be distributed equally
    % between left and right channel, ie., 50% left, 50% right, as
    % requested by the stereoPan = 0.5 setting.
    % See "PsychPortAudio DirectInputMonitoring?" for more info and details
    % on further parameters. Please note that the following commented out
    % call...
    % diResult = PsychPortAudio('DirectInputMonitoring', pa, 1);    
    % ...would do the same, the extended call is just to illustrate some
    % available optional parameters and their default settings.
    %
    % The return value 'diResult' will be zero if the call was successfull.
    % A non-zero value of diResult signals failure to enable and configure
    % direct input monitoring, either because your hardware doesn't support
    % it, or because our driver doesn't support it. See the help for
    % meaning of the different non-zero return codes.
    diResult = PsychPortAudio('DirectInputMonitoring', pa, 1, -1, 0, 0.0, 0.5);

    % Does direct input monitoring work (diResult == 0)? Then we're set.
    % Otherwise it doesn't work and we'll need to use PsychPortAudio's
    % full-duplex monitoring mode instead, which is the 2nd best
    % alternative, although certainly not zero latency:
    if diResult > 0
        % Failed! Need to use our fallback implementation:    
        fprintf('Full-duplex monitoring mode active.\n');
        PsychPortAudio('Start', pa, 0, 0, 1);
        while ~KbCheck
            WaitSecs(0.5);
            s=PsychPortAudio('GetStatus', pa);
            disp(s);
            if s.CaptureStartTime > 0
                % This Info should be available on ASIO hardware, and probably
                % Linux, but isn't available on OS/X yet:
                fprintf('Estimated minimal roundtrip latency is %f msecs.\n', 1000 * (s.StartTime - s.CaptureStartTime));
            end
        end
        PsychPortAudio('Stop', pa);
    else
        % Direct input monitoring active :-) Don't need to do anything
        % here, except waiting for a keypress from the user to finish the
        % demo:
        fprintf('Zero latency direct input monitoring mode active.\n');
        KbPressWait;
        
        % User wants us to finish. Disable input monitoring:
        PsychPortAudio('DirectInputMonitoring', pa, 0);
    end
    
    % Done - Close device and driver:
    PsychPortAudio('Close');
    
    return;
end

if ~duplex
    % Open the default audio device [], with mode 2 (== Only audio capture),
    % and a required latencyclass of latmode == low-latency mode, as well as
    % a frequency of freq Hz and 2 sound channels for stereo capture.
    % This returns a handle to the audio device:
    painput = PsychPortAudio('Open', [], 2, latmode, freq, 2, [], minLatency);
else
    % Same procedure, but open for full-duplex operation:
    painput = PsychPortAudio('Open', [], 2+1, latmode, freq, 2, [], minLatency);
    % Output- and input device are the same...
    paoutput = painput;
end

% Preallocate an internal audio recording  buffer with a capacity of at least
% 10 seconds, possibly more if requested lat'ency is higher:
PsychPortAudio('GetAudioData', painput, max(2 * lat, 10));

if ~duplex
    % Open default audio device [] for playback (mode 1), low latency (2), freq Hz,
    % stereo output:
    paoutput = PsychPortAudio('Open', [], 1, latmode, freq, 2, [], minLatency);
end

% Allocate a zero-filled (ie. silence) output audio buffer of more than
% sufficient size: Three times the requested latency, but at least 30 seconds.
% One could do this more clever, but this is a safe no-brainer and memory
% is cheap:
outbuffersize = floor(freq * 3 * max(lat, 10));
PsychPortAudio('FillBuffer', paoutput, zeros(2, outbuffersize));

% Start audio playback immediately, wait for the start to happen. Retrieve the
% start timestamp, ie., the system time when the first sample in the output
% buffer will hit the speaker in the variable 'playbackstart'.
%
% In full-duplex mode, starting the paoutput device would also start the
% painput device, as they are the same (see above), but the returned
% timestamp is always the one of start of audio output:
playbackstart = PsychPortAudio('Start', paoutput, 0, 0, 1);

% In non-duplex mode we need to start the input device separately after
% starting the output device:
if paoutput ~= painput
    % Start audio capture immediately and wait for the capture to start.
    % We set the number of 'repetitions' to zero, i.e. record/play until
    % manually stopped.
    PsychPortAudio('Start', painput, 0, 0, 1);
end

% This flag will indicate failure to achieve the wanted sound onset timing
% / latency. An experiment script would abort or reject a trial with a
% non-zero timingfailed flag:
timingfailed = 0;

% Wait until at least captureQuantum seconds of sound are available from the capture
% device and then quickly fetch it from the capture device. captureQuantum
% is the minimum amount of sound data that the driver can capture. If you'd
% ask for less you'd get at least this amount anyway + possibly extra
% delays:
s = PsychPortAudio('GetStatus', painput);
headroom = 1;
headroom = round(headroom);
captureQuantum = headroom * (s.BufferSize / s.SampleRate);
if verbose > 1
    fprintf('CaptureQuantum (Duty cycle length) is %f msecs, for a buffersize of %i samples.\n', captureQuantum * 1000, s.BufferSize);
end

[audiodata offset overflow capturestart] = PsychPortAudio('GetAudioData', painput, [], captureQuantum);

% Sanity check returned values: audiodata should be at least headroom * s.BufferSize
% samples, offset should be zero as this is the first 'GetAudioData' call
% since 'Start' of capture. overflow should be zero, otherwise we screwed
% up our timing already in the first few milliseconds because the system is
% not up to the task / overloaded for the requested latency settings.
% 'capturestart' contains the estimated time when the first returned audio
% sample hit the microphone / line-in connector:
if (size(audiodata, 2) < headroom * s.BufferSize) || (offset~=0) || (overflow > 0)
    fprintf('WARNING: SOUND ONSET TIMING SCREWED!! THE SYSTEM IS NOT UP TO THE TASK/OVERLOADED!\n');
    fprintf('Realsize samples %i < Expected size %i? Or offset %i ~= 0 ? Or overflow %i > 0 ?\n', size(audiodata, 2), headroom * s.BufferSize, offset, overflow);
    timingfailed = 1;
end

% Ok, we have our initial batch of audio samples in 'audiodata', recorded
% at time 'capturestart'. The sound output is currently feeding zeroes
% (=silence) from the zero-filled output buffer to the speakers and the
% first zero-sample in that buffer will hit the speakers at time
% 'playbackstart'. We now need to copy our 'audiodata' batch of samples
% into the output buffer, but at an offset from the start that is selected
% to exactly achieve output of our first 'audiodata' sample at the
% requested latency.
%
% The first sample was captured at time 'capturestart' and the requested
% latency for output is 'lat': Therefore the wanted playback time for this
% first sample is...
reqonsettime = capturestart + lat;

% Sanity check: Are we ahead of the playback stream with our requested
% onset time of reqonsettime? If not, then the system won't be able to
% achieve the requested 'lat'ency and we'll be late!
s = PsychPortAudio('GetStatus', paoutput);
if s.CurrentStreamTime > reqonsettime
    fprintf('WARNING: SOUND ONSET TIMING SCREWED!! THE SYSTEM IS NOT UP TO THE TASK/OVERLOADED!\n');
    fprintf(['Requested onset at time %f seconds, but audio stream is already at time %f seconds\n--> ' ...
            'We will be at least %f msecs too late!\n'], reqonsettime, s.CurrentStreamTime, 1000 * (s.CurrentStreamTime - reqonsettime));
    timingfailed = 2;
end

% The first sample from the output buffer will playback at time
% 'playbackstart', therefore our first sample should be placed at a
% timeoffset relative to the start of the outputbuffer of...
reqtimeoffset = reqonsettime - playbackstart;

% Our first audio sample needs to be placed at a time offset of
% 'reqtimeoffset' in the audio output buffer, overwriting the "silence"
% there. Map offset in seconds to offset in samples: The system plays out
% s.SampleRate samples per second, so we need to place our audio at an
% offset of...
reqsampleoffset = round(reqtimeoffset * s.SampleRate);

if reqsampleoffset < 0
    fprintf('If sound feedback works at all, then extra latency will be at least %f msecs, probably more!\n', 1000 * abs(reqtimeoffset));    
end

% Make sure the offset is positive, ie at least zero:
reqsampleoffset = max(reqsampleoffset, 0);

% Overwrite the output buffer with our captured audiodata, starting at
% sample index 'reqsampleoffset'. Need to set the 'streamingrefill' flag to
% 1 in order to enable this special overwrite mode. The 'underflow' flag
% will tell us if we made the refill in time, or if we "missed the train"
% in the last microsecond: A non-zero value means we missed.
[underflow, nextSampleStartIndex, nextSampleETASecs] = PsychPortAudio('FillBuffer', paoutput, audiodata, 1, reqsampleoffset);

s = PsychPortAudio('GetStatus', paoutput);
if underflow > 0
    fprintf('WARNING: SOUND ONSET TIMING SCREWED!! THE SYSTEM IS NOT UP TO THE TASK/OVERLOADED!\n');
    fprintf(['Requested onset at time %f seconds, but audio stream is already at time %f seconds\n--> ' ...
            'We will lose at least the first %f msecs of the sound signal!\n'], reqonsettime, s.CurrentStreamTime, 1000 * (s.CurrentStreamTime - reqonsettime));
    timingfailed = 3;
end

% Ok, if we made it until here without a non-zero 'timingfailed' flag, then
% at least the first few milliseconds of captured sound should play at
% exactly the desired 'lat'ency between capture and playback.

% From now on we'll just need to periodically fetch chunks of audio data
% from the capture device and feed it into the output device without any
% complex math or tricks involved. However in order to avoid dropouts and
% other audible artifacts we need to make sure that we feed new data fast
% enough. We will now execute a loop that tries to fetch audio in the
% smallest possible quantity from the capturedevice, then immediately
% append it to the output buffer:
updateQuantum = s.BufferSize / s.SampleRate;

% Get current status of outputdevice:
s1 = PsychPortAudio('GetStatus', paoutput);

oldcaptureQuantum = -1;
cumoverrun   = 0;
cumunderflow = 0;

% Feedback loop: Runs until ESCape keypress ...
while ~KbCheck

    % Try to dynamically adapt the amount of sound data that needs to be
    % fetched in each loop iteration. We fetch and process in larger chunks
    % if we have enough headroom. Fetching in larger 'captureQuantum'
    % chunks allows the driver to "sleep" for a few milliseconds between
    % iterations within 'GetAudioData', thereby reducing the load on the
    % operating system and cpu. This is mostly needed on MS-Windows with
    % its highly deficient scheduling and timing systems:

    % OK, this doesn't work glitch free for some reason. Leave it disabled.
    % Produces higher load, but at least works without artifacts.
    %
    % TODO: FIXME properly!
    %
    %    maxtimeUntilRefill = (nextSampleETASecs - s1.CurrentStreamTime) / 10;
    %     if maxtimeUntilRefill > (3 * updateQuantum)
    %         captureQuantum = (floor(maxtimeUntilRefill / updateQuantum) - 1) * updateQuantum;
    %     else
    %        captureQuantum = updateQuantum;
    %     end

    captureQuantum = updateQuantum;
    
    if captureQuantum ~= oldcaptureQuantum
        oldcaptureQuantum = captureQuantum;
        if verbose > 1
            fprintf('Duty cycle adapted to %f msecs...\n', 1000 * captureQuantum);
        end
    end
    
    % Get new captured sound data ...
    fetchDelay = GetSecs;
    [audiodata, offset, overrun] = PsychPortAudio('GetAudioData', painput, [], captureQuantum);
    fetchDelay = GetSecs - fetchDelay;
    
    underflow = 0;
    
    % ... and stream it into our output buffer:
    while size(audiodata, 2) > 0
        % Make sure to never push more data in the buffer than it can
        % actually hold, ie not more than half its maximum capacity:
        fetch = min(size(audiodata, 2), floor(outbuffersize / 2));
        % We feed data in chunks of 'fetch' samples:
        pushdata = audiodata(:, 1:fetch);
        % audiodata is the remainder which will be pushed in the next loop
        % iteration:
        audiodata = audiodata(:, fetch+1:end);
        
        % Perform streaming buffer refill. As long as we don't push more
        % than a buffer size, the driver will take care of the rest...
        [curunderflow, nextSampleStartIndex, nextSampleETASecs] = PsychPortAudio('FillBuffer', paoutput, pushdata, 1);
        underflow = underflow + curunderflow;
    end
    
    % Check for xrun conditions from low-level sound hardware:
    s1 = PsychPortAudio('GetStatus', paoutput);
    s2 = PsychPortAudio('GetStatus', painput);
    xruns = s1.XRuns + s2.XRuns;
    
    % Any dropouts or other audible artifacts?
    if ((overrun + underflow + xruns) > 0) && (timingfailed == 0)
        if verbose > 0
            fprintf('WARNING: SOUND DROPOUTS! THE SYSTEM IS NOT UP TO THE TASK/OVERLOADED!\n');
            fprintf('Run %i: Overruns of capture buffer: %i. Underruns of audio output buffer: %i. Hardware xruns = %i\n', tc, overrun, underflow, xruns);
        end
        timingfailed = 4;
    else
        % fprintf('nextSampleETA - currentStreamtime: %f msecs.\n', 1000 * (nextSampleETASecs - s1.CurrentStreamTime));
    end
    
    cumoverrun = cumoverrun + overrun;
    cumunderflow = cumunderflow + underflow;
    
    % Log some timing samples:
    tc = tc + 1;
    if tc <= size(tc, 2)
        tstats(:, tc) = [ s1.ElapsedOutSamples ; s1.CurrentStreamTime ; fetchDelay; nextSampleETASecs - s1.CurrentStreamTime];
    end
    
    % Done. Next iteration...
end

% Reenable all keys for KbCheck:
RestrictKeysForKbCheck([]);

% Stop the playback engine:
PsychPortAudio('Stop', paoutput, 1);

% Non-Duplex operation with separate input device?
if painput ~= paoutput
    % Stop the capture engine:
    PsychPortAudio('Stop', painput, 1);
end

% Drain its capture buffer...
PsychPortAudio('GetAudioData', painput);

% Ok, done. Close all engines and exit.
PsychPortAudio('Close');

if timingfailed > 0
    % There was trouble during execution:
    fprintf('There were timingproblems or audio dropouts during the demo [Condition %i]! Your system is not capable of reliable operation at a\n', timingfailed);
    fprintf('requested roundtrip feedback latency of %f msecs.\n\n', 1000 * lat);
    fprintf('\nOverruns of capture buffer: %i. Underruns of audio output buffer: %i. Hardware xruns = %i\n', cumoverrun, cumunderflow, xruns);
else
    fprintf('Requested roundtrip feedback latency of %f msecs seems to have worked. Please double-check with external equipment.\n\n', 1000 * lat);
end

% Prune tstats to valid range:
fprintf('Total of %i timesamples.\n', tc);
tstats = tstats(:, 1:tc);
tstats(2,:) = tstats(2,:) - tstats(2,1);
tstats(1,:) = tstats(1,:) - tstats(1,1);
[tout(1,:), idx] = unique(tstats(1,:));
tout(2:4,:) = tstats(2:4,idx);
tstats = tout;

% Plot it:
plot(tstats(1,:), tstats(2,:) * 1000, '.', tstats(1,:), tstats(3,:) * 1000, '-', tstats(1,:), tstats(4,:) * 1000, '-');

% Done.
fprintf('Demo finished, bye!\n');

return;
