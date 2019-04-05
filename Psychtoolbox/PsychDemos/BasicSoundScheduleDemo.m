function BasicSoundScheduleDemo(wavfilenames, device)
% BasicSoundScheduleDemo([wavfilenames][, device])
%
% This demo shows two things:
%
% - How to use audio schedules in PsychPortAudio to preprogram a sequence
% of different sounds to play and how to dynamically add new sounds to the
% schedule while it is playing. This is similar to "playlists" in typical
% audio player applications like iTunes or the iPod etc.
%
% - How to create and use many prefilled audio buffers before start of a
% session. This way you can preload all needed sounds before start of an
% experiment into memory, in a format optimized for fast playback and low
% memory usage. This is similar to the concept of textures or offscreen
% windows in the domain of Screen() for the visuals.
%
% Optional arguments:
%
% wavfilenames = Name of a .wav sound file to load and play back, or a cell
% array with multiple filenames to load. Otherwise the default sound files
% provided with Psychtoolbox will be loaded.
%
% device = Device index of the sound card to use.
%
%
% The demo first loads all soundfiles, and resamples them to identical
% samplingrate if possible. Then it plays the first second of each of them.
% Then it goes into an interactive mode: By pressing any of the F1 - F10
% keys or any letter key, you can select a specific file for playback. By
% pressing any other key you exit the interactive loop. After the
% interactive loop has finished, a subset of the soundfiles is played again
% with a different method, for about 2 repetitions. Then the demo exits.
%

% History:
% 04/09/2009  mk Written.
% 18-Jul-2015 mk Remove switch of 'RunMode' to 1, as 1 is the default since
%                quite a while.
% 07-Nov-2018 mk Auto-Select playback frequency, allow device selection, cosmetic.

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

% Establish mapping of keys for the interactive loop:
KbName('UnifyKeyNames');
for i=1:10
    key(i) = KbName(sprintf('F%i', i)); %#ok<AGROW>
end
for i=11:11+25
    key(i) = KbName(sprintf('%s', char('a' + i - 11))); %#ok<AGROW>
end

keyspace = KbName('space');

% Filenames provided?
if nargin < 1
    wavfilenames = [];
end

if isempty(wavfilenames)
    % No sound file provided. Load our standard sounds:
    sounddir = [ PsychtoolboxRoot 'PsychDemos' filesep 'SoundFiles' filesep ];

    % Ok, on MK's machine we have a special treat ;-)
    if exist([PsychHomeDir 'Music/StarTrekSounds/'], 'dir')
        sounddir = [ PsychHomeDir 'Music/StarTrekSounds/' ];
    end

    infilenames = dir( [ sounddir '*.wav' ] );
    infilenames.name;

    for i=1:length(infilenames)
        wavfilenames{i} = [ sounddir infilenames(i).name ];
    end
end
nfiles = length(wavfilenames);

if nargin < 2
    device = [];
end

% Always init to 2 channels, for the sake of simplicity:
nrchannels = 2;

% Perform basic initialization of the sound driver:
InitializePsychSound(1);

suggestedLatencySecs = [];
if IsARM
    % ARM processor, probably the RaspberryPi SoC. This can not quite handle the
    % low latency settings of a Intel PC, so be more lenient:
    suggestedLatencySecs = 0.025;
    fprintf('Choosing a high suggestedLatencySecs setting of 25 msecs to account for lower performing ARM SoC.\n');
end

% Open the audio 'device' with default mode [] (== Only playback),
% and a required latencyclass of 1 == standard low-latency mode, as well as
% default playback frequency and 'nrchannels' sound output channels.
% This returns a handle 'pahandle' to the audio device:
pahandle = PsychPortAudio('Open', device, [], 1, [], nrchannels, [], suggestedLatencySecs);

% Get what freq'uency we are actually using for playback:
s = PsychPortAudio('GetStatus', pahandle);
freq = s.SampleRate;

% Does a function for resampling exist? Need to load 'signal' toolbox
% on Octave for that:
if IsOctave
    try
        pkg load signal;
    catch
    end
end

if exist('resample')
    % Yes: Resample to playback sampling rate 'freq' if neccessary:
    doresample = 1
else
    % Wav files with deviating frequencies from 'freq' will play too fast or too
    % slow, because we can't resample:
    doresample = 0
end

% Read all sound files and create & fill one dynamic audiobuffer for
% each read soundfile:
buffer = [];
j = 0;

for i=1:nfiles
    try
        % Make sure we don't abort if we encounter an unreadable sound
        % file. This is achieved by the try-catch clauses...
        [audiodata, infreq] = psychwavread(char(wavfilenames(i)));
        dontskip = 1;
    catch
        fprintf('Failed to read and add file %s. Skipped.\n', char(wavfilenames(i)));
        dontskip = 0;
        psychlasterror
        psychlasterror('reset');
    end

    if dontskip
        j = j + 1;

        if doresample
            % Resampling supported. Check if needed:
            if infreq ~= freq
                % Need to resample this to target frequency 'freq':
                fprintf('Resampling from %i Hz to %i Hz... ', infreq, freq);
                audiodata = resample(audiodata, freq, infreq);
            end
        end

        [samplecount, ninchannels] = size(audiodata);
        audiodata = repmat(transpose(audiodata), nrchannels / ninchannels, 1);

        buffer(end+1) = PsychPortAudio('CreateBuffer', [], audiodata); %#ok<AGROW>
        [fpath, fname] = fileparts(char(wavfilenames(j)));
        fprintf('Filling audiobuffer handle %i with soundfile %s ...\n', buffer(j), fname);
    end
end

% Recompute number of available sounds:
nfiles = length(buffer);

% Enable use of sound schedules: We create a schedule of default size,
% currently 128 slots by default. From now on, the driver will not play
% back the sounds stored via PsychPortAudio('FillBuffer') anymore. Instead
% you'll have to define a "playlist" or schedule via subsequent calls to
% PsychPortAudio('AddToSchedule'). Then the driver will process that
% schedule by playing all defined sounds in the schedule, one after each
% other, until the end of the schedule is reached. You can add new items to
% the schedule while the schedule is already playing.
PsychPortAudio('UseSchedule', pahandle, 1);

% Build an initial play sequence. Play each buffer once for a starter:
for i=1:nfiles
    % Play buffer(i) from startSample 0.0 seconds to endSample 1.0 
    % seconds. Play one repetition of each soundbuffer...
    PsychPortAudio('AddToSchedule', pahandle, buffer(i), 1, 0.0, 1.0, 1);
end

% Suppress character spilling into command window:
ListenChar(-1);

fprintf('\nReady. Press any key to start...\n\n\n');

% Wait for single keystroke:
KbStrokeWait;

% Start audio playback of the defined schedule. We don't spec the
% 'repetitions' parameter, as this parameter is ignored when using sound
% schedules. To repeat a schedule, you must refill it in time with new
% slots or select its size so it auto-repeats properly. Why this restriction?
% As a safety measure against programming errors. Without this, it would be
% easy to create infinite playback loops due to small programming mistakes,
% where you can't exit sond playback anymore and have to kill Matlab/Octave
% just to get your system back. This restriction may be lifted in future
% driver releases if we find a good way to implement some child-protection.
%
% Start playback immediately (0) and wait for the playback to start:
PsychPortAudio('Start', pahandle, [], 0, 1);

fprintf('Audio playback started: Press key F1 to F10 or any letter key for a sound, any other key to quit.\n');

% Stay in a little interactive loop:
notprinted = 1;
while 1
    % Query current playback status:
    s = PsychPortAudio('GetStatus', pahandle);

    % Keyboard input?
    [down, secs, keyCode] = KbCheck;
    if down
        % Yes. Respond to it:
        kc = min(find(keyCode)); %#ok<MXFND>
        kc = find(key == kc);
        if ~isempty(kc)
            % A key in array 'key' pressed (F1 to F10 or a letter key).
            % Map it to a audio buffer handle:
            kc = mod(kc - 1, nfiles) + 1;
            fprintf('KEY %s pressed: Playing buffer %i\n', KbName(min(find(keyCode))), buffer(kc)); %#ok<MXFND>
            KbReleaseWait;

            % Engine still running on a schedule?
            if s.Active == 0
                % Schedule finished, engine stopped. Before adding new
                % slots we first must delete the old ones, ie., reset the
                % schedule:
                PsychPortAudio('UseSchedule', pahandle, 2);
            end

            % Add new slot with playback request for user-selected buffer
            % to a still running or stopped and reset empty schedule. This
            % time we select one repetition of the full soundbuffer:
            PsychPortAudio('AddToSchedule', pahandle, buffer(kc), 1, 0, [], 1);

            % If engine has stopped, we need to restart:
            if s.Active == 0
                PsychPortAudio('Start', pahandle, [], 0, 1);
                notprinted = 1;
            end
        else
            % Space key pressed and engine stopped?
            if keyCode(keyspace) && (s.Active == 0)
                % Reset schedule to ready state, and restart sound playback of
                % the whole current schedule:
                PsychPortAudio('UseSchedule', pahandle, 3);
                PsychPortAudio('Start', pahandle, [], 0, 1);
                notprinted = 1;
                KbReleaseWait;
            else
                % Other (unknown) key pressed: Finish interactive loop.
                break;
            end
        end
    else
        % No key pressed. Give some feedback to user:
        if (s.Active == 0) && (notprinted == 1)
            fprintf('Audio playback paused: Press key F1 to F10 or any letter key for a sound, SPACE for replay, or any other key to quit.\n');
            notprinted = 0;
        end
    end

    % Wait a bit before next status and key query. The 'YieldSecs' option
    % tells WaitSecs that this wait doesn't need to be accurate down to the
    % millisecond, but allows for some lenience in timing. This slack
    % allows to reduce system load:
    WaitSecs('YieldSecs', 0.05);
end

% End of interactive loop.
ListenChar(0);

% Stop playback: Stop immediately, but wait for stop to happen:
PsychPortAudio('Stop', pahandle, 0, 1);

% Another little demo application of audio buffers, now non-interactive:
fprintf('\n\nFillBuffer test: Using multibuffers without schedule, but as sources for streaming refill...\n');

% Disable and delete schedule: Back to standard single playbuffer operation,
% where one can only specify sounds for playback by 'FillBuffer' or
% 'RefillBuffer' commands to fill the single default playback buffer:
PsychPortAudio('UseSchedule', pahandle, 0);

% Fill playbuffer with content of buffer(1):
PsychPortAudio('FillBuffer', pahandle, buffer(1));

% Start playback in 2 seconds from now, 2 repetitions, wait for start:
PsychPortAudio('Start', pahandle, 2, GetSecs + 2, 1);

% Streaming refill with content of buffer(2). Append the content to the
% currently playing sound stream:
PsychPortAudio('FillBuffer', pahandle, buffer(2), 1);

% Streaming refill with content of buffer(4):
PsychPortAudio('FillBuffer', pahandle, buffer(4), 1);

% Wait for end of playback, then stop:
PsychPortAudio('Stop', pahandle, 1);

% Delete all dynamic audio buffers:
PsychPortAudio('DeleteBuffer');

% Close audio device, shutdown driver:
PsychPortAudio('Close');
fprintf('\n\nDone. Bye!\n\n');

return;
