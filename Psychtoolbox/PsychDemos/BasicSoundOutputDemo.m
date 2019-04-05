function BasicSoundOutputDemo(repetitions, wavfilename, device)
% BasicSoundOutputDemo([repetitions=0][, wavfilename][, device])
%
% Demonstrates very basic use of the Psychtoolbox sound output driver
% PsychPortAudio(). PsychPortAudio is a better, more reliable, more accurate
% replacement for the old Psychtoolbox SND() function and other means of
% sound output in Matlab like sound(), soundsc(), wavplay(), audioplayer()
% etc.
%
% This demo only demonstrates normal operation, not the low-latency mode.
% Other demos and tests for low-latency and high precision timing output do that.
% follow soon. If you need low-latency, make sure to read "help
% InitializePsychSound" carefully or contact the forum.
% Testing for low-latency mode showed that sub-millisecond accurate sound
% onset and < 10 msecs latency are possible on Linux, OSX and Windows with some
% sound cards.
%
%
% Optional arguments:
%
% repetitions = Number of repetitions of the sound. Zero = Repeat forever
% (until stopped by keypress), 1 = Play once, 2 = Play twice, ....
%
% wavfilename = Name of a .wav sound file to load and playback. Otherwise
% the good ol' handel.mat file (part of Matlab) is used.
%
% device = Device index of sound card to use.
%
% The demo just loads and plays the soundfile, waits for a keypress to stop
% it, then quits.

% History:
% 06/07/2007 Written (MK)
% 11/07/2018 Allow 'device' selection. (MK)

% Running on PTB-3? Abort otherwise.
AssertOpenGL;

if nargin < 1 || isempty(repetitions)
    repetitions = 0;
end

% Filename provided?
if nargin < 2
    wavfilename = [];
end

if isempty(wavfilename)
    % Ok, assign this as default sound file: Better than ol' handel - we're
    % sick of that sound.
    wavfilename = [ PsychtoolboxRoot 'PsychDemos' filesep 'SoundFiles' filesep 'funk.wav'];
end

% Read WAV file from filesystem:
[y, freq] = psychwavread(wavfilename);
wavedata = y';
nrchannels = size(wavedata,1); % Number of rows == number of channels.

% Make sure we have always 2 channels stereo output.
% Why? Because some low-end and embedded soundcards
% only support 2 channels, not 1 channel, and we want
% to be robust in our demos.
if nrchannels < 2
    wavedata = [wavedata ; wavedata];
    nrchannels = 2;
end

if nargin < 3
    device = [];
end

% Perform basic initialization of the sound driver:
InitializePsychSound;

% Open the  audio device, with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
try
    % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', device, [], 0, freq, nrchannels);
catch
    % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
    fprintf('Sound may sound a bit out of tune, ...\n\n');

    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', device, [], 0, [], nrchannels);
end

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, wavedata);

% Start audio playback for 'repetitions' repetitions of the sound data,
% start it immediately (0) and wait for the playback to start, return onset
% timestamp.
t1 = PsychPortAudio('Start', pahandle, repetitions, 0, 1);

% Wait for release of all keys on keyboard:
KbReleaseWait;

fprintf('Audio playback started, press any key for about 1 second to quit.\n');

lastSample = 0;
lastTime = t1;

% Stay in a little loop until keypress:
while ~KbCheck
    % Wait a seconds...
    WaitSecs(1);

    % Query current playback status and print it to the Matlab window:
    s = PsychPortAudio('GetStatus', pahandle);
    % tHost = GetSecs;

    % Print it:
    fprintf('\n\nAudio playback started, press any key for about 1 second to quit.\n');
    fprintf('This is some status output of PsychPortAudio:\n');
    disp(s);

    realSampleRate = (s.ElapsedOutSamples - lastSample) / (s.CurrentStreamTime - lastTime);
    fprintf('Measured average samplerate Hz: %f\n', realSampleRate);

    tHost = s.CurrentStreamTime;
    clockDelta = (s.ElapsedOutSamples / s.SampleRate) - (tHost - t1);
    clockRatio = (s.ElapsedOutSamples / s.SampleRate) / (tHost - t1);
    fprintf('Delta between audio hw clock and host clock: %f msecs. Ratio %f.\n', 1000 * clockDelta, clockRatio);
end

% Stop playback:
PsychPortAudio('Stop', pahandle);

% Close the audio device:
PsychPortAudio('Close', pahandle);

% Done.
fprintf('Demo finished, bye!\n');
