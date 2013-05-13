function BasicSoundFeedbackDemo(reqlatency, duplex)
% BasicSoundFeedbackDemo([reqlatency=7.5 ms] [,duplex=0])
%
% THIS IS EARLY ALPHA CODE! IT MAY OR MAY NOT WORK RELIABLY ON YOUR SETUP!
% TEST IT WITH MEASUREMENT EQUIPMENT IF YOU DEPEND ON ACCURATE FEEDBACK
% TIMING!!!
%
% Demonstrates very basic usage of the new Psychtoolbox sound driver
% PsychPortAudio() for audio feedback.
%
% Sound is captured from the default recording device and then - with a
% selectable delay - played back via the default output device.
%
% By default, feedback is tried with a latency of 7.5 ms plus hardware and
% system inherent delay, but you can ask for a specific latency in msecs by
% providing the optional parameter 'reqlatency'. Achievable latency will be
% constrained by the capabilities of your hardware. Choosing too low of a
% value will create audible artifacts in the sound and the driver may
% output warning about 'buffer underflows during streaming refill...'.
%
% Depending on your sound hardware you'll have to either leave 'duplex' at
% its default of zero (2 times half-duplex mode) or set it to 1
% (full-duplex mode): ASIO driven hardware -- typically on MS-Windows --
% will usually need full-duplex mode. On Macintosh OS/X it depends on the
% sound hardware. IntelMacs are happy with half-duplex mode, some PowerMacs
% may need full-duplex mode.
%
%
% If you need low-latency, make sure to read "help InitializePsychSound"
% carefully or contact the forum.
%

% History:
% 07/20/2007 Written (MK)

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

fprintf('\n\nTHIS IS EARLY ALPHA CODE! IT MAY OR MAY NOT WORK RELIABLY ON YOUR SETUP!\nTEST IT WITH MEASUREMENT EQUIPMENT IF YOU DEPEND ON ACCURATE FEEDBACK\nTIMING!!!\n\n');

% Latency provided?
if nargin < 1
    reqlatency = [];
end

if isempty(reqlatency)
    reqlatency = 7.5;
end

if nargin < 2
    duplex = [];
end

if isempty(duplex)
    duplex =0;
end

lat =  reqlatency / 1000;

% Try a sample rate of 48kHz. Should be supported by most hardware:
freq = 48000;

% Wait for release of all keys on keyboard:
KbReleaseWait;

% Perform low-level initialization of the sound driver:
InitializePsychSound(1);

% Provide some debug output:
PsychPortAudio('Verbosity', 10);

if (reqlatency == 0) && duplex
    % Special case: Full-duplex mode with minimum latency. We bypass Matlab
    % by activating PsychPortAudios full-duplex monitoring mode. The driver
    % itself will feed back all captured sound to the outputs with lowest
    % possible latency. However we don't have any control over latency or
    % sound and this only works on full-duplex hardware...
    fprintf('Full-duplex monitoring mode active.\n');
    pa = PsychPortAudio('Open', [], 4+2+1, 2, freq, 2, [], 0.010);
    PsychPortAudio('Start', pa, 0, 0, 1);
    while ~KbCheck
        WaitSecs(0.5);
        s=PsychPortAudio('GetStatus', pa);
        disp(s);
    end
    PsychPortAudio('Stop', pa);
    PsychPortAudio('Close');
    return;
end

if ~duplex
    % Open the default audio device [], with mode 2 (== Only audio capture),
    % and a required latencyclass of 2 == low-latency mode, as well as
    % a frequency of freq Hz and 2 sound channels for stereo capture.
    % This returns a handle to the audio device:
    painput = PsychPortAudio('Open', [], 2, 2, freq, 2);
else
    % Same procedure, but open for full-duplex operation:
    painput = PsychPortAudio('Open', [], 2+1, 2, freq, 2, [], 0.010);
    % Output- and input device are the same...
    paoutput = painput;
end

% Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
PsychPortAudio('GetAudioData', painput, 10);

if ~duplex
    % Open default audio device [] for playback (mode 1), low latency (2), freq Hz,
    % stereo output:
    paoutput = PsychPortAudio('Open', [], 1, 2, freq, 2);
end

% Full duplex mode. Doesn't work yet...
if duplex
    % Prefill playback buffer with silence...
    PsychPortAudio('FillBuffer', paoutput, zeros(2, freq * 2 * lat));
end

% Start audio capture immediately and wait for the capture to start.
% We set the number of 'repetitions' to zero, i.e. record/play until
% manually stopped.
painputstart = PsychPortAudio('Start', painput, 0, 0, 1);

% Wait for at least lat secs of sound data to become available: This
% directly defines a lower bound on real feedback latency. Its also a weak
% point, because waiting longer than 'lat' will increase output latency...
s=PsychPortAudio('GetStatus', painput);
availsecs = s.RecordedSecs;
while availsecs < lat
    WaitSecs(0.0001);
    s=PsychPortAudio('GetStatus', painput);
    availsecs = s.RecordedSecs;
end

% Quickly readout available sound and initialize sound output buffer with it:
[audiodata offset overflow capturestart]= PsychPortAudio('GetAudioData', painput);

if ~duplex
    % Feed everything into the initial sound output buffer:
    PsychPortAudio('FillBuffer', paoutput, audiodata);

    % Start the playback engine immediately and wait for start, let it run
    % until manually stopped:
    playbackstart = PsychPortAudio('Start', paoutput, 0, 0, 1);
else
    % Duplex mode: We don't get separate timestamps for capture and
    % playback yet. Set them to be the same - The best we can do for now...
    playbackstart = painputstart;
end

% Now the playback engine should output the first lat msecs of our sound,
% while the capture engine captures the next msecs. Compute expected
% latency. This is what the driver thinks, accuracy depends on the quality
% of implementation of the underlying sound subsystem, so its dependent on
% operating system and sound driver/sound hardware:
expecteddelay = (playbackstart - capturestart) * 1000;
fprintf('Expected latency at least %f msecs.\n', expecteddelay);

% Feedback loop: Runs until keypress ...
while ~KbCheck
    % Sleep about lat/2 secs to give the engines time to at least capture and
    % output lat/2 secs worth of sound ...
    WaitSecs(lat/5);
    
    % Get new captured sound data ...
    [audiodata, offset, overrun] = PsychPortAudio('GetAudioData', painput);
    
    % ... and stream it into our output buffer:
    while size(audiodata, 2) > 0
        % Make sure to never push more data in the buffer than it can
        % actually hold, ie not more than half its maximum capacity:
        fetch = min(size(audiodata, 2), floor(freq * lat/2));
        % We feed data in chunks...
        pushdata = audiodata(:, 1:fetch);
        audiodata = audiodata(:, fetch+1:end);
        % Perform streaming buffer refill. As long as we don't push more
        % than a buffer size, the driver will take care of the rest...
        PsychPortAudio('FillBuffer', paoutput, pushdata, 1);
    end
    
    % Done. Next iteration...
end

% Done. Stop the capture engine:
PsychPortAudio('Stop', painput, 1);

% Drain its capture buffer...
[audiodata offset]= PsychPortAudio('GetAudioData', painput);

if ~duplex
    % Stop the playback engine:
    PsychPortAudio('Stop', paoutput, 1);
end

% Ok, done. Close engines and exit.
PsychPortAudio('Close');

% Done.
fprintf('Demo finished, bye!\n');

return;
