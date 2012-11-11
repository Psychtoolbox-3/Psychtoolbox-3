function varargout = DatapixxAudioKey(cmd, varargin)
% DatapixxAudioKey - A prototype implementation of a DataPixx
% audio-key/voice-key.
%
% This allows to record audio from the Datapixx, optionally synchronized to
% visual stimulus onset, optionally with simple but precise timestamping of
% onset of audio signals.
%
% This is an experimental proof-of-concept prototype! Its name, behaviour
% and interface may change without warning in future beta-releases!
%
%
% Subfunctions and their meaning:
% -------------------------------
%
% DatapixxAudioKey('Open' [, sampleRate][, lrMode][, inputJack][, inputGain]);
% - Open audio system for recording. This will stop any running audio
% acquisition operations. Audio sampling will be performed at the given
% optional 'sampleRate' between 8000 Hz and 96Khz with the given number of
% audio input channels 'lrMode' (0 = Mono: Average of left and right
% channel, 1 = Only left channel, 2 = Only right channel, 3 = Stereo).
% Optionally select 'inputJack' as audio input (1=Microphone, 2=Line-In)
% and 'inputGain' as input amplifier gain. If no 'inputJack' is specified
% and no input is preselected, 'inputJack' will default to microphone
% input.
%  
%
% DatapixxAudioKey('Close');
% - Close audio subsystem after stopping any running acquisition operations.
%
% [startTimeBox, startTimeGetSecs] = DatapixxAudioKey('CaptureNow' [, maxDurationSecs=30][, startOffsetSecs=0][,bufferBaseAddress=20e6] [, numBufferFrames=maxScheduleFrames]);
% - Start audio capture immediately (ie., with minimum possible delay on
% your system), return a 'startTimeBox' timestamp in Datapixx clock time of
% when capture actually started or will start, taking the 'startOffsetSecs'
% into account. You can also use the 2nd optional return argument
% 'startTimeGetSecs' to get the same timestamp mapped to Psychtoolbox
% GetSecs() time. For a more precise mapping to GetSecs time, you can
% convert the 'startTimeBox' timestamp into Psychtoolbox GetSecs() time
% with the proper mapping functions of PsychDataPixx (see "help
% PsychDataPixx").
%
% The optional 'startOffsetSecs' allows to delay the start of the capture
% operation by 'startOffsetSecs' seconds. The optional 'maxDurationSecs'
% allows to define an upper limit onto the duration of the audio capture
% operation. The operation will stop automatically after the given number
% of seconds. By default, capture will run for 30 seconds. For
% robustness and efficiency, you should specify a reasonable upper limit if
% possible. 'bufferBaseAddress' is the memory start address of the audio
% buffer for capture, 'numBufferFrames' is the size of the buffer. See
% "Datapixx SetAudioSchedule?" for more info. Leave at default settings if
% possible.
%
%
% DatapixxAudioKey('CaptureAtFlip'[, flipCount=next][, maxDurationSecs=30][, startOffsetSecs=0][,bufferBaseAddress=20e6] [, numBufferFrames=maxScheduleFrames]);
% - Schedule start of audio capture synchronized to the visual stimulus
% onset of a future Screen('Flip') or Screen('AsyncFlipBegin') command.
% All parameters are identical to the ones for DatapixxAudioKey('CaptureNow',...),
% except for the first optional parameter 'flipCount'. 'flipCount' defines
% at which Screen('Flip') audio capture should be started. If omitted or
% set to zero or [], capture will be triggered by the next flip command.
% Otherwise it will be triggered by the flip command with the given
% 'flipCount'. You can query the current flipCount via PsychDatapixx('FlipCount').
%
% Example: Setting flipCount = PsychDatapixx('FlipCount') + 10; would start
% audio capture at the 10'th invocation of a flip command from now.
%
% This function doesn't return a capture start timestamp as capture will
% only happen in the future, so the timestamp is only available in the
% future.
%
%
% [audiodata, onsetTimeSecs] = DatapixxAudioKey('GetResponse'[, amountSecs=current][, blocking=1][, stopCapture=0]);
% - Try to fetch available captured audiodata from a running capture
% operation. Only call this function after a capture operation has been
% started via DatapixxAudioKey('CaptureNow') or scheduled for start at a
% certain flipCount via DatapixxAudioKey('CaptureAtFlip') and the flip is
% imminent, i.e., after the flip command for that flip has been called
% already. The function will error-out if you call it too early!
%
% The optional 'amountSecs' asks the driver to return exactly
% 'amountSecs' worth of audio data. By default it will return whatever
% is available at time of call. If 'amountSecs' is specified, but the
% requested amount is not yet available, the optional 'blocking' flag will
% define behaviour: If set to 1 (default), the driver will wait until the
% specified amount becomes available. If set to 0, the driver will return
% with empty [] return arguments so you can retry later. The optional
% 'stopCapture' flag if set to 1 will stop audio capture, its default
% setting of zero will keep capture running until manually stopped or until
% it stops by itself.
%
% 'audiodata' is the vector or matrix of captured audiodata. 1 Row for mono
% recording, or a 2-row matrix (one row for each audio channel) in stereo
% recording modes. Each value is an audio signal sample in range [-1 ; 1].
%
% 'onsetTimeSecs' is the time since start of audio capture when a certain
% level of loudness (as set by DatapixxAudioKey('TriggerLevel')) was first
% exceeded, e.g., due to onset of a vocal response of a subject.
%
% 'readOffset' is the sample index of the next sample that will be read on
% the next call to this function, if any.
%
%
% oldLevel = DatapixxAudioKey('TriggerLevel' [, newLevel]);
% - Return old and optionally set new trigger threshold level for the
% timestamping of onset of audio signals in DatapixxAudioKey('GetResponse').
%
% 'oldLevel' is the current/old level. 'newLevel' is the optional new
% level. Level can be between 0 and 1, with a default level of 0.1 for 10%
% of max signal intensity as trigger level.
%
%
% DatapixxAudioKey('AutoTriggerLevel', silenceSecs);
% - Auto-Select trigger threshold level for audio onset timestamping. After
% start of audio capture, there must be a period of silence, where the only
% signal recorded is due to noise in the microphone/amplifiers etc, not due
% to actual signal. You must provide the minimum duration of this time
% period of silence as 'silenceSecs' parameter, e.g., 1.5 for 1.5 seconds
% of guaranteed silence at start of recording. During recording, the
% function will use this chunk of silence to compute an optimal trigger
% level, which is 10% louder than the loudest sample in the silence block
% and then use that triggerlevel as if DatapixxAudioKey('TriggerLevel') had
% been called with that auto-selected level.
%
%
% DatapixxAudioKey('StopCapture');
% - Stop audio capture as soon as possible.
%
%

persistent dpixaudioin;

if nargin < 1 || isempty(cmd)
    error('Required Subcommand missing or empty!');
end

if isempty(dpixaudioin) || (dpixaudioin.refcount == 0)
    % Driver closed at startup:
    dpixaudioin.refcount = 0;
    dpixaudioin.samplerate = 0;
    dpixaudioin.pendingForFlip = -1;
    dpixaudioin.triggerLevel = 0.1;
    dpixaudioin.readOffset = 0;
    dpixaudioin.lrMode = 0;
end

if strcmpi(cmd, 'Open')
    % Open the audiokey: Performs basic init and setup of audio input
    % subsystem.
    if dpixaudioin.refcount > 0
        error('Open: Tried to open DatapixxAudioKey, but audio key already open!');
    end
    
    % Open device connection if not already open:
    PsychDataPixx('Open');
    
    if PsychDataPixx('IsBusy')
        error('Open: Tried to open DatapixxAudioKey, but Datapixx is busy! Screen flip pending?');
    end
    
    % Stop potentially running audio input schedule:
    Datapixx('StopMicrophoneSchedule');
    
    % TODO FIXME: How to decide if this is needed or harmful?
    Datapixx('InitAudio');
    
    % Disable loopback:
    Datapixx('DisableAudioLoopback');
    
    % Execute:
    Datapixx('RegWrRd');
    
    micstatus = Datapixx('GetMicrophoneStatus');

    % Assign new optional sampleRate:
    if length(varargin) >= 1 && ~isempty(varargin{1})
        dpixaudioin.samplerate = varargin{1};
    else
        dpixaudioin.samplerate = micstatus.scheduleRate;
    end

    % Assign new optional lrMode:
    if length(varargin) >= 2 && ~isempty(varargin{2})
        dpixaudioin.lrMode = varargin{2};
    else
        dpixaudioin.lrMode = micstatus.lrMode;
    end

    % Assign new optional input source:
    if length(varargin) >= 3 && ~isempty(varargin{3})
        inputjack = varargin{3};
    else
        % Leave current input setting "as is":
        inputjack = micstatus.source;

        % If no input set at all, default to microphone input:
        if inputjack == 0
            inputjack = 1;
        end
    end
    
    % Assign new optional input gain:
    if length(varargin) >= 4 && ~isempty(varargin{4})
        inputgain = varargin{4};
    else
        inputgain = micstatus.gain;
    end

    % Assign audio input source and input gain:
    Datapixx('SetMicrophoneSource', inputjack, inputgain);
    Datapixx('RegWrRd');

    % Mark us as open:
    dpixaudioin.refcount = dpixaudioin.refcount + 1;
    
    return;    
end

if strcmpi(cmd, 'Close')
    if dpixaudioin.refcount <= 0
        % Noop:
        dpixaudioin.refcount = 0;
        return;
    end
    
    if PsychDataPixx('IsBusy')
        error('Close: Tried to close DatapixxAudioKey, but Datapixx is busy! Screen flip pending?');
    end
    
    % Stop potentially running audio input schedule:
    Datapixx('StopMicrophoneSchedule');
    
    % Execute:
    Datapixx('RegWrRd');
    
    % Close our connection:
    PsychDataPixx('Close');
    
    % Closed:
    dpixaudioin.refcount = 0;
    
    return;
end

% All following commands need the audiokey to be open:
if dpixaudioin.refcount <= 0
    error('Tried to use audiokey, but audiokey not open!');
end

if strcmpi(cmd, 'CaptureNow')
    if dpixaudioin.pendingForFlip > -1
        error('CaptureNow: Tried to start or engage audiokey, but audiokey already started or engaged!');
    end

    if length(varargin) >= 1 && ~isempty(varargin{1})
        maxScheduleFrames = varargin{1} * dpixaudioin.samplerate;
    else
        % Default to 30 seconds if nothing provided:
        maxScheduleFrames = 30 * dpixaudioin.samplerate;
    end

    if length(varargin) >= 2 && ~isempty(varargin{2})
        scheduleOffsetSecs = varargin{2};
    else
        scheduleOffsetSecs = 0;
    end

    if PsychDataPixx('IsBusy')
        error('CaptureNow: Datapixx is busy! Screen flip pending?');
    end
    
    % Stop potentially running audio input schedule:
    Datapixx('StopMicrophoneSchedule');
    Datapixx('RegWrRd');
    
    % Init audio input schedule:
    if length(varargin) >= 3
        Datapixx('SetMicrophoneSchedule', scheduleOffsetSecs, dpixaudioin.samplerate, maxScheduleFrames, dpixaudioin.lrMode, varargin{3:end});
    else
        Datapixx('SetMicrophoneSchedule', scheduleOffsetSecs, dpixaudioin.samplerate, maxScheduleFrames, dpixaudioin.lrMode);
    end

    % Reset readOffset in audio capture buffer:
    dpixaudioin.readOffset = 0;
    
    % Start audio input schedule:
    Datapixx('StartMicrophoneSchedule');

    % Execute:
    Datapixx('RegWrRd');
    
    if nargout > 0
        % Return boxtime of start of capture:
        tStartBox = Datapixx('GetTime');
        varargout{1} = tStartBox + scheduleOffsetSecs;
    end
    
    if nargout > 1
        % Return GetSecs mapped time as well:
        tStartGetSecs = PsychDataPixx('FastBoxsecsToGetsecs', tStartBox);
        varargout{2} = tStartGetSecs;
    end
    
    % Mark capture as started:
    dpixaudioin.pendingForFlip = 0;

    return;
end

if strcmpi(cmd, 'CaptureAtFlip')
    if dpixaudioin.pendingForFlip > -1
        error('CaptureAtFlip: Tried to start or engage audiokey, but audiokey already started or engaged!');
    end
    
    if length(varargin) < 1 || isempty(varargin{1}) || varargin{1} == 0
        targetFlip = [];
    else
        targetFlip = varargin{1};
    end
    
    if targetFlip <= PsychDataPixx('FlipCount')
        % targetFlip set to a count that is already over. Pointless! Set it
        % to "next flip" instead:
        targetFlip = [];
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        maxScheduleFrames = varargin{2} * dpixaudioin.samplerate;
    else
        % Default to 30 seconds if nothing provided:
        maxScheduleFrames = 30 * dpixaudioin.samplerate;
    end

    if length(varargin) >= 3 && ~isempty(varargin{3})
        scheduleOffsetSecs = varargin{3};
    else
        scheduleOffsetSecs = 0;
    end
    
    if PsychDataPixx('IsBusy')
        error('CaptureAtFlip: Datapixx is busy! Screen flip pending?');
    end
    
    % TODO FIXME: Should stop and configuration of schedule be part of the
    % command packet to be executed at flip? Would provide the ability to
    % "preprogram" complex sequences of audio response collection, but make
    % overall design much more complex?
    %
    % For now we go with the simple solution:

    % Stop potentially running audio input schedule:
    Datapixx('StopMicrophoneSchedule');
    
    % Init audio input schedule:
    if length(varargin) >= 4
        Datapixx('SetMicrophoneSchedule', scheduleOffsetSecs, dpixaudioin.samplerate, maxScheduleFrames, dpixaudioin.lrMode, varargin{4:end});
    else
        Datapixx('SetMicrophoneSchedule', scheduleOffsetSecs, dpixaudioin.samplerate, maxScheduleFrames, dpixaudioin.lrMode);
    end

    % Execute:
    Datapixx('RegWrRd');

    % Reset readOffset in audio capture buffer:
    dpixaudioin.readOffset = 0;
    
    % Build command string for psync-triggered start of schedule at target
    % flipcount:
    cstr = 'Datapixx(''StartMicrophoneSchedule'');';
    
    % Schedule cstr for execution at stimulus onset for the 'targetFlip'th
    % flip command execution:
    PsychDataPixx('ExecuteAtFlipCount', targetFlip, cstr);
    dpixaudioin.pendingForFlip = targetFlip;
    
    % Ready.
    return;
end

if strcmpi(cmd, 'StopCapture')
    if PsychDataPixx('IsBusy')
        error('StopCapture: Datapixx is busy! Screen flip pending?');
    end
    
    % Submit stop command immediately:
    Datapixx('StopMicrophoneSchedule');
    Datapixx('RegWrRd');
    
    % Reset flag:
    dpixaudioin.pendingForFlip = -1;
    
    return;
end

if strcmpi(cmd, 'TriggerLevel')
    varargout{1} = dpixaudioin.triggerLevel;
    if length(varargin) >= 1 && ~isempty(varargin{1})
        dpixaudioin.triggerLevel = varargin{1};
    end
    return;
end

if strcmpi(cmd, 'AutoTriggerLevel')
    % Assign value in seconds as negative number:
    dpixaudioin.triggerLevel = -1 * varargin{1};
    return;
end

if strcmpi(cmd, 'GetResponse')
    if dpixaudioin.pendingForFlip < 0
        error('GetResponse: Tried to get response from audiokey, but audiokey not active or engaged!');
    end

    % If audio recording was scheduled to start in sync with a certain
    % Screen flip, check if that specific target flip count has been
    % reached:
    if (dpixaudioin.pendingForFlip > 0) && (PsychDataPixx('FlipCount') < dpixaudioin.pendingForFlip)
        % Audio key shall start recording at a certain flipcount which
        % likely hasn't been reached yet (according to PsychDataPixx('FlipCount')).
        %
        % Multiple options:
        %
        % 1 We're more than one flip away from startpoint. In that case no
        %   point waiting for start, as usercode would first need to
        %   execute some Screen('Flip/AsyncFlipBegin') command to make it
        %   even possible to start audio capture.
        %
        % 2 We're exactly one flip away and...
        %   a) No flip command submitted by usercode -> No way this is
        %   gonna fly.
        %
        %   b) Asyncflip command submitted: Flip will eventually happen and
        %   get us going. This is the case if the device is marked as "busy
        %   waiting for psync" and async flip is marked as active.
        %
        % In case 1 or 2a we can't wait for a response, whereas in case 2b
        % we can wait for a response as it will eventually happen.
        
        % Check for case 1:
        if PsychDataPixx('FlipCount') < dpixaudioin.pendingForFlip - 1
            % More than 1 flip away. That's a no-no:
            fprintf('DatapixxAudioKey: Start of audio input scheduled for flipcount %i, but flipcount not yet reached.\n', dpixaudioin.pendingForFlip);
            fprintf('DatapixxAudioKey: Current count %i is more than 1 flip away! Impossible to "GetResponse" from device this way!\n', PsychDataPixx('FlipCount'));
            error('GetResponse: Tried to get response from audiokey, but audiokey not yet active and impossible to activate by this command-sequence!');            
        end
        
        % Check for case 2b:
        if PsychDataPixx('IsBusy') && Screen('GetWindowInfo', PsychDataPixx('WindowHandle'), 4)
            % Asyncflip with Psync op pending. The flip will eventually
            % happen at some point in the future, so we can enter a polling
            % loop to wait for it to happen.
        else
            % 1 flip away but no flip scheduled to actually reach that
            % trigger-flip:
            fprintf('DatapixxAudioKey: Start of audio input scheduled for flipcount %i, but flipcount not yet reached!\n', dpixaudioin.pendingForFlip);
            fprintf('DatapixxAudioKey: Current count %i. Impossible to "GetResponse" from device this way,\n', PsychDataPixx('FlipCount'));
            fprintf('DatapixxAudioKey: because no flip has been scheduled to reach the trigger flipcount!\n');
            error('GetResponse: Tried to get response from audiokey, but audiokey not yet active and impossible to activate by this command-sequence!');
        end
    end

    if length(varargin) < 2 || isempty(varargin{2})
        blocking = 1;
    else
        blocking = varargin{2};
    end
    
    % Ok, should be already running or it is at least certain that it will
    % eventually run due to trigger by psync mechanism.
    %
    % Enter a polling loop to wait until the requested minimum amount of
    % audio data is available:
    if length(varargin) >= 1 && ~isempty(varargin{1})
        % Wait or poll for requested amount of audiodata:
        numFrames = ceil(varargin{1} * dpixaudioin.samplerate);
        while 1
            Datapixx('RegWrRd');
            micstatus = Datapixx('GetMicrophoneStatus');
            if micstatus.newBufferFrames >= numFrames
                % Requested amount of audio data available: Exit loop.
                break;
            else
                % Insufficient amount. If this is a polling request, we
                % simply return no result:
                if ~blocking
                    varargout{1} = [];
                    varargout{2} = [];
                    return;
                end
            end
            
            % Sleep a msec, then retry:
            WaitSecs('YieldSecs', 0.001);
        end
    else
        % Just get current status - current amount of audiodata:
        Datapixx('RegWrRd');
        micstatus = Datapixx('GetMicrophoneStatus');
        numFrames = -1;
    end
    
    % Either requested amount of audio data available or no specific amount
    % requested. Stop schedule if requested by usercode:
    if length(varargin) >= 3 && ~isempty(varargin{3}) && varargin{3} == 1
        Datapixx('StopMicrophoneSchedule');
        Datapixx('RegWrRd');
        dpixaudioin.pendingForFlip = -1;
        micstatus = Datapixx('GetMicrophoneStatus');
    end

    % Fetch all captured frames if no specific amount requested:
    if numFrames == -1
        numFrames = micstatus.newBufferFrames;
    end

    % Retrieve Datapixx internal delay, so we can correct for it:
    inDelaySecs = Datapixx('GetMicrophoneGroupDelay', dpixaudioin.samplerate);

    % Retrieve recorded data:
    audiodata = Datapixx('ReadMicrophoneBuffer', ceil(numFrames), -1);
    
    % Compute timestamp in seconds since start of capture of when the
    % triggerLevel was exceeded the first time:
    if dpixaudioin.triggerLevel >= 0
        % Specific triggerLevel given. Use "as is":
        triggerTime = (dpixaudioin.readOffset + min(find(abs(audiodata(1,:)) > dpixaudioin.triggerLevel))) / dpixaudioin.samplerate; %#ok<MXFND>
    else
        % Auto-Trigger mode. Leadtime of silence in seconds given. We can
        % use that many seconds of signal at the beginning, knowing it
        % represents silence, to determine an optimal triggerLevel. We set
        % a level that is 10% more than the "loudest" signal sample in the
        % "silence" lead-in block of samples:
        triggerLevel = 1.1 * max(abs(audiodata(1, 1:floor(dpixaudioin.samplerate * abs(dpixaudioin.triggerLevel)))));
        
        if triggerLevel > 0.9
            warning('Psychtoolbox:DatapixxAudioKey:AutoTriggerHigh', 'Auto-Selected triggerLevel %f exceeds 90% of max possible signal amplitude. Strong noise or something wrong with your setup?!?', triggerLevel);
        end
        
        % Do trigger threshold seeking like above:
        triggerTime = (dpixaudioin.readOffset + min(find(abs(audiodata(1,:)) > triggerLevel))) / dpixaudioin.samplerate; %#ok<MXFND>        
    end
    
    % Compensate for Datapix hw delays:
    triggerTime = triggerTime - inDelaySecs;
    
    % Increment readOffset into audio capture buffer:
    dpixaudioin.readOffset = dpixaudioin.readOffset + size(audiodata, 2);
    
    % Return all data:
    varargout{1} = audiodata;
    varargout{2} = triggerTime;
    varargout{3} = dpixaudioin.readOffset;
    
    return;
end

error('DatapixxAudioKey: Unknown subcommand provided!');
