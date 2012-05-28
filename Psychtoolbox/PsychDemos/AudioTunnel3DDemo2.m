function AudioTunnel3DDemo2(idx, fdir)
% AudioTunnel3DDemo2 -- A very sketchy and raw demo of how to use OpenAL to
% create the illusion of walking through a tunnel which has different sound
% sources attached to it.
%
% You'll need a properly configured multi-speaker system or very good
% head-phones attached to a very good soundcard and a lot of luck and
% goodwill for this to sound good.
%
% OpenAL for Matlab is currently only supported on a subset of platforms,
% so this demo may abort quickly with an error.
%

if nargin < 1
    idx = 1;
end

if nargin < 2
    fdir = 1;
end

startpos = -20;
walkspeed = 12;
maxbehind = 4;
radius = 2;

nsources = 1;

% Establish key mapping: ESCape aborts, Space toggles between auto-
% movement of sound source or user mouse controlled movement:
KbName('UnifyKeynames');
space = KbName('space');
esc = KbName('ESCAPE');

% Initialize OpenAL subsystem at debuglevel 2 with the default output device:
InitializeMatlabOpenAL(2);

% Generate one sound buffer:
buffers = alGenBuffers(nsources);

% Query for errors:
alGetString(alGetError)

% Try to load some impressive sound...
sounddir = [PsychtoolboxRoot 'PsychDemos/SoundFiles/'];
soundfiles = dir([sounddir '*.wav'])

alListenerfv(AL.POSITION, [0, 0, 0]);
alListenerfv(AL.VELOCITY, [0, 0, 0]);

if IsOSX
    alcASASetListener(ALC.ASA_REVERB_ON, 1);
    alcASASetListener(ALC.ASA_REVERB_QUALITY, ALC.ASA_REVERB_QUALITY_Max);
    alcASASetListener(ALC.ASA_REVERB_ROOM_TYPE, ALC.ASA_REVERB_ROOM_TYPE_Cathedral);
end

% Create a sound source:
sources = alGenSources(nsources);

perm = randperm(nsources);

for i=1:nsources
    % Assign next soundfilename round-robin:
    soundname = [sounddir soundfiles(mod(perm(i), length(soundfiles))+1).name];
    soundname = [sounddir soundfiles(mod(idx + i, length(soundfiles))+1).name];
    
    % Load it...
    [mynoise freq]= wavread(soundname);
    mynoise = mynoise(:, 1);
    
    % Convert it...
    mynoise = int16(mynoise * 32767);

    % Fill our sound buffer with the data from the sound vector. Tell AL that its
    % a 16 bpc, mono format, with length(mynoise)*2 bytes total, to be played at
    % a sampling rate of freq Hz. The AL will resample this to the native device
    % sampling rate and format at buffer load time.
    alBufferData( buffers(i), AL.FORMAT_MONO16, mynoise, length(mynoise)*2, freq);

    % Attach our buffer to it: The source will play the buffers sound data.
    alSourceQueueBuffers(sources(i), 1, buffers(i));

    % Switch source to looping playback: It will repeat playing the buffer until
    % its stopped.
    alSourcei(sources(i), AL.LOOPING, AL.TRUE);

    % Set emission volume to 100%, aka a gain of 1.0:
    alSourcef(sources(i), AL.GAIN, 1);

%     alSourcef(sources(i), AL.CONE_INNER_ANGLE, 30);
%     alSourcef(sources(i), AL.CONE_OUTER_ANGLE, 270);
    alSource3f(sources(i), AL.DIRECTION, 0, 0, 1);

    % Set position: Nicely lined up in z direction, but randomly placed in
    % x,y direction within a range of +/- 5 meters:
    z(i) = maxbehind + 1;
    
    % Sources themselves remain static in space:
    alSource3f(sources(i), AL.VELOCITY, 0, 0, fdir * walkspeed);
    
    if IsOSX
        % Source emits some sound that gets reverbrated in room:
        alcASASetSource(ALC.ASA_REVERB_SEND_LEVEL, sources(i), 0.0);
    end
end

% Start playback for these sources:
alSourcePlayv(nsources, sources);

% 3D sound animation loop. Runs until ESCape key press:
manual = 0;
tstart = GetSecs;

while 1
    % Check keyboard:
    [isdown dummy, keycode]=KbCheck;
    if isdown
        if keycode(esc)
            break;
        end

        if keycode(space)
            % Toggle between mouse control and automatic:
            manual=1-manual;
            KbReleaseWait;
        end
    end
    
    if manual
        % Query mouse:
        [xm ym]=GetMouse;
        walkspeed=(ym-500)/100;
    end

    for i=1:nsources
        if z(i) > maxbehind
            theta = rand * 2 * pi;
            x(i) = cos(theta) * radius;
            y(i) = sin(theta) * radius;
            z(i) = startpos;
        end

        alSource3f(sources(i), AL.POSITION, x(i), y(i), fdir * z(i));
    end

    t = GetSecs;
    telapsed = t - tstart;
    tstart = t;
    tdistance = telapsed * walkspeed;
    z = z + tdistance;

    % Pause for 10 milliseconds in order to yield the cpu to other processes:
    WaitSecs(0.01);
end

% Stop playback of all sources:
alSourceStopv(length(sources), sources);

for i=1:nsources
    % Unqueue sound buffer:
    alSourceUnqueueBuffers(sources(i), 1, buffers(i));
end

% Wait a bit:
WaitSecs(0.1);

% Delete buffer:
alDeleteBuffers(nsources, buffers);

% Wait a bit:
WaitSecs(0.1);

% Delete sources:
alDeleteSources(nsources, sources);

% Wait a bit:
WaitSecs(0.1);

% Shutdown OpenAL:
CloseOpenAL;

% Done. Bye.
return;
