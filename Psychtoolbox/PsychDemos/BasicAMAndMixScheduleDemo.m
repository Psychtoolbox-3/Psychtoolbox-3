function BasicAMAndMixScheduleDemo
% Demonstrates basic use of sound schedules, volume controls, amplitude
% modulation and audio mixing of multiple voices.
%
% Usage: BasicAMAndMixScheduleDemo;
%
% For more advanced or different use of sound schedules see
% BasicSoundScheduleDemo. For more advanced use of mixing, modulators and
% other stuff see BasicMultiTrackAudioRecorderDemo.
%
% This demo shows how to use master-slave audio devices to playback with
% multiple voices, mixing multiple sounds together. It shows how to
% manually control the per-slave volume and master volume. It shows how to
% apply AM amplitude modulation to the different slave audio devices, and how
% to synchronize different devices and AM modulation. It also shows use of
% sound schedules to implement repetitive sound sequences with precise
% timing between different sounds.
%

% History:
% 25.04.2010   mk  Written.


KbName('UnifyKeyNames');
leftArrow = KbName('LeftArrow');
rightArrow = KbName('RightArrow');
upArrow = KbName('UpArrow');
downArrow = KbName('DownArrow');
escape = KbName('ESCAPE');

% Initialize Sounddriver:
InitializePsychSound(1);

% Open real default [] soundcard as master device (+8) for playback only (+1), with
% standard low-latency, high timing precision mode, 2 channels, 48kHz:
nrchannels = 2;
freq = 48000;

% Add 15 msecs latency on Windows, to protect against shoddy drivers:
sugLat = [];
if IsWin
    sugLat = 0.015;
end

pamaster = PsychPortAudio('Open', [], 1+8, 1, freq, nrchannels, [], sugLat);

% Start master immediately, wait for it to be started. We won't stop the
% master until the end of the session.
PsychPortAudio('Start', pamaster, 0, 0, 1);

% Set the masterVolume for the master: This volume setting affects all
% attached sound devices. We set this to 0.5, so it doesn't blow out the
% ears of our listeners...
PsychPortAudio('Volume', pamaster, 0.5);

% Create to slave audio devices for sound playback (+1), with same
% frequency, channel count et. as master. Attach them to master. As they're
% attached to the same sound channels of the master (actually the same
% single channel), their audio output will mix together:
pasound1 = PsychPortAudio('OpenSlave', pamaster, 1);
pasound2 = PsychPortAudio('OpenSlave', pamaster, 1);

% Attach some beep-tone to each of them. pasound1 gets a 500 Hz, 1 sec
% tone, pasound2 gets a 750 Hz, 1 sec tone:
sound1 = MakeBeep(500, 1.0, freq);
sound2 = MakeBeep(750, 1.0, freq);
sound3 = MakeBeep(1000, 1.0, freq);

% Make sure we have always 2 channels stereo output.
% Why? Because some low-end and embedded soundcards
% only support 2 channels, not 1 channel, and we want
% to be robust in our demos.
sound1 = [sound1 ; sound1];
sound2 = [sound2 ; sound2];
sound3 = [sound3 ; sound3];

% Create audio buffers prefilled with the 3 sounds:
pabuffer1 = PsychPortAudio('CreateBuffer', [], sound1);
pabuffer2 = PsychPortAudio('CreateBuffer', [], sound2);
pabuffer3 = PsychPortAudio('CreateBuffer', [], sound3);

% Demo of sound schedules: Build a sound schedule (a "playlist") that will
% play the three sound buffers infinitely repeating. Each sound is played 2
% seconds after start of playback of the previous one. As the sounds
% themselves are 1 second in duration, you should hear 1st sound for 1 sec,
% then 1 sec silence, then 2nd sound for a second, 1 sec silence, 3rd
% sound for a second, 1 sec silence, 1st sound ...
% The schedule shall contain exactly 6 slots which will repeat ad
% infinitum, or until stopped:
PsychPortAudio('UseSchedule', pasound1, 1, 6);

% This command code in a slot tells to take a break (+1) before processing
% the following of the following slot. The (+16) means to wait until the
% given amount of seconds has elapsed since start of playback of the most
% recent soundbuffer. Therefore it defines a relative spacing between
% playback of successive sound buffers.
% Note that there are more command codes available (see PsychPortAudio AddToSchedule?).
% These allow for other types of timing, just don't make sense in this
% demo.
cmdCode = 1 + 16;

% Add 1st sound buffer: The special '1' flag at the end tells not to delete
% this command from the schedule, but keep it for repeated execution on
% repetitions of the schedule:
PsychPortAudio('AddToSchedule', pasound1, pabuffer1, [], [], [], [], 1);

% Tell pasound1 to start playing the following buffer exactly 2.0 seconds
% after playback of the previous buffer has started:
PsychPortAudio('AddToSchedule', pasound1, -cmdCode, 2.0, [], [], [], 1);

% Add 2nd sound buffer:
PsychPortAudio('AddToSchedule', pasound1, pabuffer2, [], [], [], [], 1);

% Tell pasound1 to start playing the following buffer exactly 2.0 seconds
% after playback of the previous buffer has started:
PsychPortAudio('AddToSchedule', pasound1, -cmdCode, 2.0, [], [], [], 1);

% Add 3rd sound buffer:
PsychPortAudio('AddToSchedule', pasound1, pabuffer3, [], [], [], [], 1);

% Tell pasound1 to start playing the following buffer exactly 2.0 seconds
% after playback of the previous buffer has started:
PsychPortAudio('AddToSchedule', pasound1, -cmdCode, 2.0, [], [], [], 1);

fprintf('Starting sound playback 3 seconds from now...\n');
tStart = PsychPortAudio('Start', pasound1, 0, GetSecs + 3, 1);
fprintf('Started! Will add some AM modulation of first tone at start of 3rd repetition...\n');

% Create an AM modulator (+32) to gain-modulate the signal of pasound1, that is
% perform amplitude modulation. Modulation will only happen while
% pamodulator1 is active. No modulation will happen if it is stopped:
pamodulator1 = PsychPortAudio('OpenSlave', pasound1, 32);

% 'envelope1' sound is 1 second at 8 Hz frequency,
% built for a playback device with 'freq' Hz sampling
% rate:
envelope1 = (1 + MakeBeep(8, 1, freq)) / 2;

% Fill it into standard sound buffer of pamodulator1, instead of using
% 'CreateBuffer', just for a change...
PsychPortAudio('FillBuffer', pamodulator1, [envelope1; envelope1]);

% We define a repeating schedule for the pamodulator1 as well. This time
% only 2 slots. The first plays the 1 second AM modulation signal in our
% 'FillBuffer'ed buffer 0, then the 2nd slot pauses until 6 seconds after
% start of our modulation signal. This way we will get a 1 second modulation
% repeating every 6 seconds. Because the period of 6 seconds of
% pamodulator1 is identical to the 6 seconds period of the sound sequence
% of pasound1, they are phase-locked, so this will basically modulate the
% same 1 second piece of sound in pasound1:
PsychPortAudio('UseSchedule', pamodulator1, 1, 2);
PsychPortAudio('AddToSchedule', pamodulator1, 0, [], [], [], [], 1);
PsychPortAudio('AddToSchedule', pamodulator1, -cmdCode, 6.0, [], [], [], 1);

% Modulation schedule ready. We want to start modulation at the 3rd overall
% repetition of the three-tone sequence that pasound1 is playing. We know
% that one three-tone sequence takes 6 seconds and that the first iteration
% of the sequence started at 'tStart', so the 3rd repetition must therefore
% begin at tThirdRep = tStart + 2 * 6  seconds. Schedule start of our
% pamodulator1 for target time tThirdRep. This should cause modulation of
% the 1st tone, starting at the 3rd overall repetition, then repeating
% until stopped:
tThirdRep = tStart + 2 * 6;

% Schedule start of AM on infinite repeat. Wait for the actual start to happen:
PsychPortAudio('Start', pamodulator1, 0, tThirdRep, 1);

fprintf('Press any key to continue.\n');
KbStrokeWait;

% Stop the modulation immediately. This will just output the sequence at
% full constant volume again:
PsychPortAudio('Stop', pamodulator1);
fprintf('Back to full-blast unmodulated sequence!\n');

% Ok, this time we create a nice sin'e curve AM modulation to "soften" the
% onset and offset of each of the three tones - nicer to your ears :-)
%
% Each AM envelope is again 1 second duration and we choose a repeating
% pattern of 2 second period. This way for each 6 seconds playout of the
% 3-tone sequence, we will get 3 repetitions of the 2 seconds modulation
% sequence, phase-locked to the 3 seconds sequence, nicely modulating each
% of the three tones:

% Redefine our envelope, softly raising from 0.0 - 1.0, then falling from
% 1.0 to 0.0 again, all in 1 second:
envelope2 = sin((0:1/freq:1) * pi);

% Fill it into standard sound buffer of pamodulator1, instead of using
% 'CreateBuffer', just for a change...
PsychPortAudio('FillBuffer', pamodulator1, [envelope2; envelope2]);

% Clear and reset our 2-slot schedule:
PsychPortAudio('UseSchedule', pamodulator1, 2);

% Our 1 second envelope:
PsychPortAudio('AddToSchedule', pamodulator1, 0, [], [], [], [], 1);

% Restart playback every 2 seconds for a total period of 2 seconds:
PsychPortAudio('AddToSchedule', pamodulator1, -cmdCode, 2.0, [], [], [], 1);

% Ok, we need to start the new schedule in sync with start of a new 3-tone
% sequence for best results. We don't know where we are in the sequence,
% but we can query when the sequence started repeating the last time. It is
% stored in the status.StartTime field of pasound1 device.
% We need to wait until playout of the first tone started again to get the
% reference timestamp we want. Poll for this event:
while 1
    status = PsychPortAudio('GetStatus', pasound1);
    % Every 6th slot corresponds again to the 1st tone, and it needs to be
    % playing that is status.Active:
    if (mod(status.SchedulePosition, 6) == 0) && (status.Active)
        % Got it! The status.startTime now corresponds to the start of the
        % current 3-tone-sequence:
        break;
    end
    
    % Nope. Wait a bit, then retry:
    WaitSecs('YieldSecs', 0.1);    
end

% Let's start our modulator offset by 2 repetitions, that is finish the
% current repetition of the 6 second cycle, then play another one, then
% start modulating. We do this here for simplicity, so we don't need to
% worry being already very close to the start of the next repetition by
% pure chance and thereby possibly missing the deadline:
tModStart = status.StartTime + 2 * 6;

PsychPortAudio('Start', pamodulator1, 0, tModStart);

fprintf('Will start soft modulation of each tone %f seconds from now at %f seconds.\n', tModStart - GetSecs, tModStart);
fprintf('This is %f seconds since start.\n', tModStart - tStart);

fprintf('Press any key to continue.\n');
KbStrokeWait;

% Ok, let's add the second pasound device into the mix. This one will play
% white noise, again envelope modulated at a low frequency:

% 5 Second worth of random noise in -1 ; +1 range:
noise = 2 * rand(1, 5 * freq) - 1;

% Into pasound2's standard buffer it goes...
PsychPortAudio('FillBuffer', pasound2, [noise ; noise]);

% Create a dedicated AM modulator for pasound2 as well:
pamodulator2 = PsychPortAudio('OpenSlave', pasound2, 32);

% Fill it with a slow version of the soft envelope of pamodulator1:
% Redefine our envelope, softly raising from 0.2 - 1.0, then falling from
% 1.0 to 0.2 again, all in 6 seconds:
envelope3 = 0.2 + 0.8 * sin((0:1/(freq*6):1) * pi);

% Fill it into standard sound buffer of pamodulator2, instead of using
% 'CreateBuffer', just for a change...
PsychPortAudio('FillBuffer', pamodulator2, [envelope3 ; envelope3]);

% Ok, we don't use any schedules for pasound2 and pamodulator2, we had that
% already. Instead we just start both of them, set to infinite repeat. We
% need to start them at exactly the same time so the modulation starts
% softly at the very first played out white-noise sample. We could pick any
% point in time in the future, e.g., GetSecs + 0.5 seconds, but for the fun
% of it, we will start it again synchronized to the sequence of pasound1,
% recycling the same trick we used before. Start it two repetitions after
% start of the most recent repetition of the whole 6-slot schedule.
% We need to wait until playout of the first tone started again to get the
% reference timestamp we want. Poll for this event:
while 1
    status = PsychPortAudio('GetStatus', pasound1);
    % Every 6th slot corresponds again to the 1st tone, and it needs to be
    % playing that is status.Active:
    if (mod(status.SchedulePosition, 6) == 0) && (status.Active)
        % Got it! The status.startTime now corresponds to the start of the
        % current 3-tone-sequence:
        break;
    end
    
    % Nope. Wait a bit, then retry:
    WaitSecs('YieldSecs', 0.1);    
end
tModStart = status.StartTime + 2 * 6;

fprintf('Will start additional modulated random noise %f seconds from now.\n', tModStart - GetSecs);
fprintf('Should start at onset of first tone of the 3-tone sequence.\n');

% Schedule start of pamodulator2 at tModStart, set it to infinite repeat:
PsychPortAudio('Start', pamodulator2, 0, tModStart);

% Schedule start of pasound2 at tModStart, set it to infinite repeat.
% We wait until start has happened and retrieve its onset time, just
% because we can...
tNoiseStart = PsychPortAudio('Start', pasound2, 0, tModStart, 1);

fprintf('Started at system time %f seconds. Feel the breeze?\n\n', tNoiseStart);
fprintf('Press ESCAPE key to continue, arrow keys to manually change relative volumes.\n');
while 1
    [secs, keyCode] = KbStrokeWait;
    if keyCode(escape)
        break;
    end

    if keyCode(leftArrow)
        vs = max(PsychPortAudio('Volume', pasound1) - 0.1, 0);
        PsychPortAudio('Volume', pasound1, vs);
    end
    
    if keyCode(rightArrow)
        vs = min(PsychPortAudio('Volume', pasound1) + 0.1, 10);
        PsychPortAudio('Volume', pasound1, vs);
    end

    if keyCode(downArrow)
        vs = max(PsychPortAudio('Volume', pasound2) - 0.1, 0);
        PsychPortAudio('Volume', pasound2, vs);
    end
    
    if keyCode(upArrow)
        vs = min(PsychPortAudio('Volume', pasound2) + 0.1, 10);
        PsychPortAudio('Volume', pasound2, vs);
    end

end

% Let's stop pasound1 exactly at the end of one of its sound cycles. Find
% start of next cycle:
while 1
    status = PsychPortAudio('GetStatus', pasound1);
    % Every 6th slot corresponds again to the 1st tone, and it needs to be
    % playing that is status.Active:
    if (mod(status.SchedulePosition, 6) == 0) && (status.Active)
        % Got it! The status.startTime now corresponds to the start of the
        % current 3-tone-sequence:
        break;
    end
    
    % Nope. Wait a bit, then retry:
    WaitSecs('YieldSecs', 0.1);    
end

% Compute stop time exactly 6 seconds after start:
tStop = status.StartTime + 6;

fprintf('Will stop tones at end of this 3-tone sequence.\n');

% Request (flag 3) a stop at time 'tStop'. Wait (flag 1) for stop to
% happen. Please note that a 'waitForEndOfPlayback' flag of 1 or use of the
% 'repetition' parameter wouldn't work here, as we've got a schedule
% attached and these parameters don't work with schedules...
PsychPortAudio('Stop', pasound1, 3, 1, [], tStop);

% We can immediately stop the modulator, as there ain't any sound anymore to
% modulate anyway:
PsychPortAudio('Stop', pamodulator1);

fprintf('3-Tone sequence over. Now softly fading out the breeze over about 10 seconds...\n');

% Fade out the breeze noise manually via pulling down the volume in 1000
% steps, 1 step every 0.01 seconds:
vs = PsychPortAudio('Volume', pasound2);
for v = vs:-0.001:0
    % Assign new reduced masterVolume to pasound2:
    PsychPortAudio('Volume', pasound2, v);
    WaitSecs(0.01);
end

% Ok, we're fully muted. Let's immediately stop both the sound device and
% its modulator:
PsychPortAudio('Stop', pasound2);
PsychPortAudio('Stop', pamodulator2);

% That's it. Everything stopped and silent. Close all devices, release all
% ressources, shutdown the driver:
PsychPortAudio('Close');

fprintf('Finished. Bye!\n');

return;

