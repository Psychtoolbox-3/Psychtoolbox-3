function BasicSoundInputDemo(wavfilename, voicetrigger, maxsecs, device)
% BasicSoundInputDemo([wavfilename][, voicetrigger=0][, maxsecs=inf] [, device])
%
% Demonstrates very basic usage of the new Psychtoolbox sound driver
% PsychPortAudio() for audio capture / recording.
%
% Sound is captured from the default recording device, optionally waiting
% until the amplitude exceeds some threshold level before start of
% recording. At the end of the recording session the recorded sound is
% played back vie the default output device. Waveform data of captured
% sound is also plotted to a Matlab figure window during capture.
%
% This demo only demonstrates normal operation, not the low-latency mode,
% extra demos and tests for low-latency and high precision timing output will
% follow soon. If you need low-latency, make sure to read "help
% InitializePsychSound" carefully or contact the forum.
%
% Optional arguments:
%
% wavfilename = Name of a .wav sound file to store the recorded sound to.
%               If left out, sound won't be stored to filesystem.
% 
% voicetrigger = If set to a non-zero threshold value, the driver will wait
%                for the sound signal to exceed the specified voicetrigger threshold
%                level before it starts capturing audio data.
%
% maxsecs      = Maximum number of seconds of sound to capture. Defaults to
%                infinite - sound is recorded until a key is pressed.
%
% device       = Deviceindex of audio card to use. Auto-Selected if omitted.

% History:
% 06/30/2007 Written (MK)
% 08/10/2008 Add some sound onset time calculation for the fun of it. (MK)
% 11/07/2018 Auto select samplerate, use low-latency, workaround Octave plotting,
%            allow 'device' selection. Cosmetic. (MK)

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

% Filename provided?
if nargin < 1
    wavfilename = [];
end

if nargin < 2
    voicetrigger = [];
end

if isempty(voicetrigger)
    voicetrigger = 0;
end

if nargin < 3
    maxsecs = [];
end

if isempty(maxsecs)
    maxsecs = inf;
end

if nargin < 4
    device = [];
end

% Workaround broken qt plotting on some Octave setups:
if IsOctave && exist('graphics_toolkit')
    try
        graphics_toolkit ('fltk');
    catch
    end
end

% Wait for release of all keys on keyboard:
KbReleaseWait;

% Perform basic initialization of the sound driver:
InitializePsychSound;

% Open audio device 'device', with mode 2 (== Only audio capture),
% and a required latencyclass of 1 == low-latency mode, with the preferred
% default sampling frequency of the audio device, and 2 sound channels
% for stereo capture. This returns a handle to the audio device:
pahandle = PsychPortAudio('Open', device, 2, 1, [], 2);

% Get what freq'uency we are actually using:
s = PsychPortAudio('GetStatus', pahandle);
freq = s.SampleRate;

% Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
PsychPortAudio('GetAudioData', pahandle, 10);

% Start audio capture immediately and wait for the capture to start.
% We set the number of 'repetitions' to zero,
% i.e. record until recording is manually stopped.
PsychPortAudio('Start', pahandle, 0, 0, 1);

fprintf('Audio capture started, press any key for about 1 second to quit.\n');

% Want to start via VoiceTrigger?
if voicetrigger > 0
    % Yes. Fetch audio data and check against threshold:
    level = 0;

    % Repeat as long as below trigger-threshold:
    while level < voicetrigger
        % Fetch current audiodata:
        [audiodata offset overflow tCaptureStart] = PsychPortAudio('GetAudioData', pahandle);

        % Compute maximum signal amplitude in this chunk of data:
        if ~isempty(audiodata)
            level = max(abs(audiodata(1,:)));
        else
            level = 0;
        end

        % Below trigger-threshold?
        if level < voicetrigger
            % Wait before next scan:
            WaitSecs(0.0001);
        end
    end

    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(audiodata(1,:)) >= voicetrigger)); %#ok<MXFND>

    % Initialize our recordedaudio vector with captured data starting from
    % triggersample:
    recordedaudio = audiodata(:, idx:end);

    % For the fun of it, calculate signal onset time in the GetSecs time:
    tOnset = tCaptureStart + ((offset + idx - 1) / freq);

    fprintf('Estimated signal onset time is %f secs, this is %f msecs after start of capture.\n', tOnset, (tOnset - tCaptureStart)*1000);
else
    % Start with empty sound vector:
    recordedaudio = [];
end

% We retrieve status once to get access to SampleRate:
s = PsychPortAudio('GetStatus', pahandle);

% Stay in a little loop until keypress:
while ~KbCheck && ((length(recordedaudio) / s.SampleRate) < maxsecs)
    % Wait a second...
    WaitSecs(1);

    % Query current capture status and print it to the Matlab window:
    s = PsychPortAudio('GetStatus', pahandle);

    % Print it:
    fprintf('\n\nAudio capture started, press any key for about 1 second to quit.\n');
    fprintf('This is some status output of PsychPortAudio:\n');
    disp(s);

    % Retrieve pending audio data from the drivers internal ringbuffer:
    audiodata = PsychPortAudio('GetAudioData', pahandle);
    nrsamples = size(audiodata, 2);

    % Plot it, just for the fun of it:
    plot(1:nrsamples, audiodata(1,:), 'r', 1:nrsamples, audiodata(2,:), 'b');
    drawnow;

    % And attach it to our full sound vector:
    recordedaudio = [recordedaudio audiodata]; %#ok<AGROW>
end

% Stop capture:
PsychPortAudio('Stop', pahandle);

% Perform a last fetch operation to get all remaining data from the capture engine:
audiodata = PsychPortAudio('GetAudioData', pahandle);

% Attach it to our full sound vector:
recordedaudio = [recordedaudio audiodata];

% Close the audio device:
PsychPortAudio('Close', pahandle);

% Replay recorded data: Open 'device' for output, push recorded sound
% data into its output buffer:
pahandle = PsychPortAudio('Open', device, 1, 0, freq, 2);
PsychPortAudio('FillBuffer', pahandle, recordedaudio);

% Start playback immediately, wait for start, play once:
PsychPortAudio('Start', pahandle, 1, 0, 1);

% Wait for end of playback, then stop engine:
PsychPortAudio('Stop', pahandle, 1);

% Close the audio device:
PsychPortAudio('Close', pahandle);

% Shall we store recorded sound to wavfile?
if ~isempty(wavfilename)
    psychwavwrite(transpose(recordedaudio), freq, 16, wavfilename)
end

% Done.
fprintf('Demo finished, bye!\n');

