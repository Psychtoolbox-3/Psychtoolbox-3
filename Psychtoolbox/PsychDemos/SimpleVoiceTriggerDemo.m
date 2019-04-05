function SimpleVoiceTriggerDemo(triggerlevel, device)
% SimpleVoiceTriggerDemo([triggerlevel=0.1][, device])
%
% Demonstrates very basic usage of the new Psychtoolbox sound driver
% PsychPortAudio() for implementation of a "voice trigger". This really
% only collects the "voice response time" it doesn't actually store the
% response itself. Have a look at BasicSoundInputDemo for how one can
% actually retrieve and save the sound vector with the response itself.
%
% In any case you *must* verify correct timing of your sound hardware with
% some external measurement equipment, e.g., in conjunction with the
% PsychPortAudioTimingTest script (or the AudioFeedbackLatencyTest for a
% crude test). You only need to do this once after major changes to your
% systems software, hardware or Psychtoolbox installation, but there isn't
% any sane way to avoid such a validation!
%
% There are two different methods demonstrated to do so. Please read the
% source code of this file for details.
%
% In method I:
%
% The "response window" is fixed to 5 seconds. A trial will always last 5
% seconds, regardless if you respond early, late, or not at all.
%
%
% In method II:
%
% Sound is captured from the default recording device, waiting
% until the amplitude exceeds some 'triggerlevel'.
%
% If the triggerlevel is exceeded, sound capture stops, returning the
% estimated time of "voice onset" in system time.
%

% History:
% 08/10/2008 Written (MK)
% 11/07/2018 Auto select samplerate, allow 'device' selection. Cosmetic. (MK)

if nargin < 1 || isempty(triggerlevel)
    triggerlevel = 0.1;
    fprintf('No "triggerlevel" argument in range 0.0 to 1.0 provided: Will use default of 0.1...\n\n');
end

if nargin < 2
    device = [];
end

% Perform basic initialization of the sound driver, initialize for
% low-latency, high timing precision mode:
InitializePsychSound(1);

% Open the audio device, with mode 2 (== Only audio capture),
% and a required latencyclass of two 2 == low-latency mode, as well as
% default frequency and 2 sound channels for stereo capture. We also
% set the required latency to a pretty high 20 msecs. Why? Because we don't
% actually need low-latency here, we only need low-latency mode of
% operation so we get maximum timing precision -- Therefore we request
% low-latency mode, but loosen our requirement to 20 msecs.
%
% This returns a handle to the audio device:
pahandle = PsychPortAudio('Open', device, 2, 2, [], 2, [], 0.02);

% Get what freq'uency we are actually using:
s = PsychPortAudio('GetStatus', pahandle);
freq = s.SampleRate;

fprintf('\n\nPART I - Simple "offline" response collection at end of trial.\n');
fprintf('Each trials response period lasts 5 seconds...\n\n');

% FIRST METHOD: Supersimple
% This is suitable if the duration of the response collection period / the
% duration of a trial is fixed and if the script does not need to respond
% in any way to voice onset -- if the "voice trigger" is actually not
% triggering anything, but just the onset timestamp is of interest.
%
% In this case we simply start sound recording at the beginning of the
% trial - into a sound buffer of fixed and sufficient capacity.
% At the end of a trial (or a response period), we stop recording, fetch
% all audio data from the buffer (as well as a timestamp of when exactly
% recording was started), then use some simple Matlab thresholding code to
% find the first sample in the sound data vector which is above threshold,
% ie., the sample of response onset. Then we translate that samples index
% into a relative time to start of recording, add that to the absolute
% start time of recording, and get the absolute system time of voice onset,
% which we can log, or compare to other timestamps from Screen('Flip'),
% GetSecs, whatever to compute reaction times.

% Preallocate an internal audio recording  buffer with a generous capacity
% of 10 seconds: We only need about 5 seconds, but let's be generous...
PsychPortAudio('GetAudioData', pahandle, 10);

% Do five trials:
for trial = 1:5
    % Start audio capture immediately and wait for the capture to start.
    % We set the number of 'repetitions' to zero, i.e. record until recording
    % is manually stopped.
    PsychPortAudio('Start', pahandle, 0, 0, 1);

    % Tell user to shout:
    fprintf('Make noise! Make noise!          ');
    tStim = GetSecs;

    % Wait for about 5 seconds, so user has time to shout:
    WaitSecs(5);

    % Stop sound capture: End of response period.
    PsychPortAudio('Stop', pahandle);

    % Fetch all about 5 seconds of audiodata at once:
    [audiodata offset overflow tCaptureStart]= PsychPortAudio('GetAudioData', pahandle);

    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(sum(audiodata)) >= triggerlevel)); %#ok<MXFND>

    % Any response?
    if isempty(idx)
        fprintf('No response at all within 5 seconds?!?\n\n');
    else
        % Compute absolute event time:
        tOnset = tCaptureStart + ((offset + idx - 1) / freq);

        % Print RT:
        fprintf('---> Reaction time is %f milliseconds.\n', (tOnset - tStim)*1000);
    end

    % Next trial after 2 seconds:
    WaitSecs(2);
end

% Done with method I. You'd now close the audio device, but we'll do that
% at the end of this script instead...

% =========================================================================

% SECOND METHOD: Slightly more elaborate
% The principle is the same as in method I, but now we don't analyze the
% whole chunk of audio data at the end of trial, but we implement some
% polling loop which periodically fetches small consecutive chunks of sound
% data and does the analysis. This way, response collection / stimulus
% presentation / whatever can be stopped as soon as we have our voice onset
% timestamp --> The voice trigger can actually really trigger some action
% in our script! This also allows for "infinite" response collection -- as
% we're only chewing on small chunks of audio data, the total duration of
% response collection is not in any way limited by available memory or
% such.
fprintf('\n\nPART II - The saga continues... This time with a polling method.\n\n');

% Do ten trials:
for trial = 1:10
    % Preallocate an internal audio recording  buffer with a generous capacity
    % of 10 seconds:
    PsychPortAudio('GetAudioData', pahandle, 10);

    % Start audio capture immediately and wait for the capture to start.
    % We set the number of 'repetitions' to zero, i.e. record until recording
    % is manually stopped.
    PsychPortAudio('Start', pahandle, 0, 0, 1);

    % Tell user to shout:
    fprintf('Make noise! Make noise!          ');
    tStim = GetSecs;

    % Wait in a polling loop until some sound event of sufficient loudness
    % is captured:
    level = 0;

    % Repeat as long as below trigger-threshold:
    while level < triggerlevel
        % Fetch current audiodata:
        [audiodata offset overflow tCaptureStart]= PsychPortAudio('GetAudioData', pahandle);

        % Compute maximum signal amplitude in this chunk of data:
        if ~isempty(audiodata)
            level = max(abs(sum(audiodata)));
        else
            level = 0;
        end

        % Below trigger-threshold?
        if level < triggerlevel
            % Wait for five milliseconds before next scan:
            WaitSecs(0.005);
        end
    end

    % Ok, last fetched chunk was above threshold!
    % Find exact location of first above threshold sample.
    idx = min(find(abs(sum(audiodata)) >= triggerlevel)); %#ok<MXFND>

    % Compute absolute event time:
    tOnset = tCaptureStart + ((offset + idx - 1) / freq);

    % Stop sound capture:
    PsychPortAudio('Stop', pahandle);

    % Fetch all remaining audio data out of the buffer - Needs to be empty
    % before next trial:
    PsychPortAudio('GetAudioData', pahandle);

    % Print RT:
    fprintf('---> Reaction time is %f milliseconds.\n', (tOnset - tStim)*1000);

    % Next trial after 2 seconds:
    WaitSecs(2);
end

% Close the audio device:
PsychPortAudio('Close', pahandle);

% Done.
fprintf('Demo finished, bye!\n');

return;
