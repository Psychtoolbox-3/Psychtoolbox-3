function SimpleSoundScheduleDemo
% SimpleSoundScheduleDemo - Show simple application of sound schedules.
%
% This demo shows a very basic application of PsychPortAudio's sound
% schedule functions. It defines a sequence of 1 second beep tones which are
% supposed to be emitted at specific times relative to some trigger event,
% e.g., 5 seconds, 10 seconds and 15 seconds after the trigger event. The
% trigger event here is a key press, but could be something else like a
% visual stimulus onset or reception of some other external trigger signal.
%
% The demo performs 3 trials, each waiting for a keypress, then executing
% the 3 beep tone schedule. After those 3 trials it exits.
%
% For more complex applications and use of sound schedules see
% BasicSoundScheduleDemo. For a complex mixture with audio mixing, see
% BasicAMAndMixScheduleDemo.
%

% History:
% 7-Aug-2013  mk  Written.

% Select playback sampling rate of 44100 Hz:
freq = 44100;

% Select one channel mono playback:
nrchannels = 1;

% Open sound device 'pahandle' with specified freq'ency and number of audio
% channels for playback in timing precision mode on the default audio
% device:
InitializePsychSound(1);
pahandle = PsychPortAudio('Open', [], [], [], freq, nrchannels);

% Fill default audio buffer with a 500 Hz beep sound of 1 second duration:
snddata = MakeBeep(500, 1, freq);
PsychPortAudio('FillBuffer', pahandle, snddata);

% Enable use of sound schedules to define "playlists" of sounds to play:
PsychPortAudio('UseSchedule', pahandle, 1);

% Define relative sound onset times, relative to 'Start' of playback:
% Play first sound 5 seconds after 'triggerTime' below, 2nd sound 10
% seconds and 3rd sound 15 seconds after 'triggerTime':
soundTimes = [5, 10, 15];
soundTimes = sort(soundTimes);

% Convert into relative delays, relative to 'triggerTime' for first sound
% soundTimes(1) and then relative to each predecessor for following sounds
% via the diff() op:
delaysound = [soundTimes(1) , diff(soundTimes)];

% Define use of timed delays, relative to requested start of previous sound:
delayCmd = -(1+8);

% Play sound in audio buffer delaySound(1) seconds after base time:
PsychPortAudio('AddToSchedule', pahandle, delayCmd, delaysound(1));
PsychPortAudio('AddToSchedule', pahandle, 0);

% Play next sound delaySound(2) seconds after start of previous sound:
PsychPortAudio('AddToSchedule', pahandle, delayCmd, delaysound(2));
PsychPortAudio('AddToSchedule', pahandle, 0);

% Play next sound delaySound(3) seconds after start of previous sound:
PsychPortAudio('AddToSchedule', pahandle, delayCmd, delaysound(3));
PsychPortAudio('AddToSchedule', pahandle, 0);

% Repeat for 3 trials:
for trial = 1:3
    % Wait for some trigger signal, could be also Screen('flip',...) for
    % visual onset, define its occurence as baseline time 'triggerTime':
    fprintf('Hit a key to start audio schedule for trial %i.\n', trial);
    triggerTime = KbStrokeWait;
    
    % Start audio schedule, with times relative to 'triggerTime':
    PsychPortAudio('Start', pahandle, [], triggerTime);
    
    % Wait until sound schedule ends playing:
    fprintf('Waiting for audio schedule to finish playing.\n');
    PsychPortAudio('Stop', pahandle, 3);
    
    % Rewind schedule so it can be used in another trial:
    PsychPortAudio('UseSchedule', pahandle, 3);
end

% End of session, close down driver:
PsychPortAudio('Close', pahandle);

return;
