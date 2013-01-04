function varargout = ResponsePixx(cmd, varargin)
% ResponsePixx - Control and use the ResponsePixx response button box.
%
% This allows to record button events from the ResponsePixx response box
% device. Start and stop of logging of button responses and control of the
% button lights (on/off and intensity) can be done immediately or
% synchronized to the operation of the Screen('Flip') command, ie.,
% synchronized to visual stimulus onset.
%
% By default, ResponsePixx response boxes with 5 buttons are supported.
% Specify the 'nrButtons' parameter in the 'Open' subfunction if you need a
% different button count.
%
% The code can scan the first 16 digital input pins of the Datapixx, i.e.,
% nrButtons can be up to 16, e.g., to control self-made response boxes. The
% last 8 input pins are not useable for timestamped button queries and are
% instead used to drive the button illumination of the ResponsePixx
% devices.
%
% Caution: Avoid intermixing 'StartNow'/'StopNow' commands with
% 'StartAtFlip'/'StopAtFlip' commands. While this may work, there is some
% potential for using it the wrong way and producing very "funny" bugs.
%
% For reaction time measurements wrt. visual stimuli, see the timestamp
% logging functions of the PsychDataPixx driver to get timestamps in
% Datapixx clock time. You can also calculate RT's in GetSecs time wrt. the
% host clock by using the usual Screen('Flip') or PsychPortAudio stimulus
% onset timestamps, which may be even more convenient. However, if you do
% so, make sure that you read and understand the comments about achieving
% good precision by correcting for clock drift, both in the help of the
% ResponsePixx subfunctions and in the help of PsychDataPixx. Specifically
% look at PsychDataPixx('GetPreciseTime') and
% PsychDataPixx('BoxsecsToGetsecs').
%
%
% Button mapping on the ResponsePixx handheld device:
% ---------------------------------------------------
%
% Buttonvector element [1,2,3,4,5] == [Red, Yellow, Green, Blue, White].
%
%
% Subfunctions and their meaning:
% -------------------------------
%
% ResponsePixx('Open' [, numSamples=1000][, bufferBaseAddress=12e6][, nrButtons=5]);
% - Open ResponsePixx for button response collection. Configures the 24 bit
% digital input port of the DataPixx to receive button responses from a
% ResponsePixx response button box and to drive the button illumination
% lights of that box. Enabled button debouncing on the digital inputs.
%
% Optional argument 'numSamples' defines the maximum number of button
% transitions to log. Defaults to 1000 if omitted. 'bufferBasAddress' if
% provided is the memory start address of the logging buffer inside the
% Datapixx. Defaults to 12e6 if omitted.
%
% Optional argument 'nrButtons' defines the number of buttons to query. Can
% be up to 16 buttons for the first 16 TTL inputs of the DataPixx, but
% defaults to 5 buttons - the maximum number supported by the current
% ResponsePixx devices.
%
% This function doesn't start logging of responses. Call the 'StartNow' or
% 'StartAtFlip' subfunction for start of response collection. Call
% 'StopNow' or 'StopAtFlip' to stop response collection.
%
%
% ResponsePixx('Close');
% - Close ResponsePixx after stopping any running logging operations.
% Disable the button lights and reset the digital input port of the
% Datapixx.
%
%
% [startTimeBox, startTimeGetSecs] = ResponsePixx('StartNow' [, clearLog=0][, buttonLightState][, buttonLightIntensity]);
% - Start button response collection immediately (ie., with minimum
% possible delay on your system). Set the optional 'clearLog' flag to 1 if
% all currently stored button responses should be discarded prior to
% logging of new responses, e.g., at start of a new trial.
%
% The optional vector 'buttonLightState' allows to define the on/off state
% of the button illumination, e.g., buttonLightState = [0,1,0,0,1] would
% turn on the light inside button 2 and 5 and turn off the lights in
% buttons 1, 3 and 4. By default, the state of the button lights is not
% changed. The optional value 'buttonLightIntensity' in the range 0.0 to
% 1.0 controls the intensity of the button illumination.
%
% Returns a 'startTimeBox' timestamp in Datapixx clock time of when
% acquisition actually started. You can also use the 2nd optional return
% argument 'startTimeGetSecs' to get the same timestamp mapped to
% Psychtoolbox GetSecs() time. For a more precise post-hoc mapping to
% GetSecs time, you can convert the 'startTimeBox' timestamp into
% Psychtoolbox GetSecs() time with the proper mapping functions of
% PsychDataPixx (see "help PsychDataPixx").
%
% The function can be called multiple times, while logging is already
% started, if you just want to clear the log of stored responses or change
% the state of the button lights.
%
%
% ResponsePixx('StartAtFlip' [, flipCount=next][, clearLog=0][, buttonLightState][, buttonLightIntensity]);
% - Schedule start of response collection synchronized to the visual stimulus
% onset of a future Screen('Flip') or Screen('AsyncFlipBegin') command.
%
% All parameters are identical to the ones for ResponsePixx('StartNow',...),
% except for the first optional parameter 'flipCount'. 'flipCount' defines
% at which Screen('Flip') acquisition should be started. If omitted or
% set to zero or [], acquisition will be started by the next flip command.
% Otherwise it will be started by the flip command with the given
% 'flipCount'. You can query the current flipCount via PsychDatapixx('FlipCount').
%
% Example: Setting flipCount = PsychDatapixx('FlipCount') + 10; would start
% acquisition at the 10'th invocation of a flip command from now.
%
% This function doesn't return a start timestamp as acquisition will
% only happen in the future, so the timestamp is only available in the
% future.
%
%
% [stopTimeBox, stopTimeGetSecs] = ResponsePixx('StopNow' [, clearLog=0][, buttonLightState][, buttonLightIntensity]);
% - Stop button response collection immediately (ie., with minimum
% possible delay on your system). See 'StartNow' for meaning of optional
% arguments.
%
% The function can be called multiple times, while logging is already
% stopped, if you just want to clear the log of stored responses or change
% the state of the button lights.
%
%
% ResponsePixx('StopAtFlip' [, flipCount=next][, clearLog=0][, buttonLightState][, buttonLightIntensity]);
% - Schedule stop of response collection synchronized to the visual stimulus
% onset of a future Screen('Flip') or Screen('AsyncFlipBegin') command.
% See 'StartAtFlip' for meaning of parameters.
%
%
% [buttonStates, transitionTimesSecs, underflows] = ResponsePixx('GetLoggedResponses' [, numberResponses=current][, blocking=1][, timeout=inf]);
% - Try to fetch logged button responses from a running logging
% operation. Only call this function after logging has been
% started via ResponsePixx('StartNow') or scheduled for start at a
% certain flipCount via ResponsePixx('StartAtFlip') and the flip is
% imminent, i.e., after the flip command for that flip has been called
% already. The function will error-out if you call it too early!
%
% The optional 'numberResponses' asks the driver to return exactly
% 'numberResponses' button transitions. By default it will return whatever
% is available at time of call. If 'numberResponses' is specified, but the
% requested amount is not yet available, the optional 'blocking' flag will
% define behaviour: If set to 1 (default), the driver will wait until the
% specified amount becomes available. If set to 0, the driver will return
% with empty [] return arguments so you can retry later. If 'blocking' is
% set, then the driver will wait for at most 'timeout' seconds for a
% response, then return anyway. By default the timeout is infinite.
%
% 'buttonStates' is a n-by-5 matrix: The n'th row encodes the button state
% after the n'th transition of button state. Each column encodes up/down
% state of a button, ie., the j'th column encodes state of the j'th button.
% 0 = Button released, 1 = Button pressed. A button press followed by a
% release would create two rows of button transition, one for the
% buttonState after pressing a button, a second one for the state after
% releasing the button. The ResponsePixx uses a debouncing algorithm to
% filter out multiple transitions within a 30 msecs time window.
%
% 'transitionTimesSecs' is the vector of timestamps in Datapixx clock time
% of when the corresponding button state has changed. The i'th element
% corresponds to the time when the i'th element in the 'buttonStates'
% vector was logged due to a change in button state at that time. See
% "help PsychDataPixx" on how to map these box timestamps to GetSecs()
% timestamps if you wish.
%
% 'underflows' is the number of times the logging buffer underflowed. This
% indicated lost or corrupted button events. It should only happen if your
% subject pressed more buttons within two calls to 'GetLoggedResponses'
% than the event buffer can hold. The size of the buffer can be selected at
% 'Open' time with the 'numSamples' parameter, but it defaults to a
% generous 1000 button events.
%
%
% [buttonStates, boxTimeSecs, getsecsTimeSecs] = ResponsePixx('GetButtons');
% - Perform immediate query of response box button states. Return current
% state in 'buttonStates' and corresponding box time and GetSecs time
% timestamps of time of query in 'boxTimeSecs' and 'getsecsTimeSecs'.
%
%

% History:
% 5.9.2010   mk  Written.
% 7.9.2010   mk  Refined, now really useable.

persistent rpixstatus;

if nargin < 1 || isempty(cmd)
    error('Required Subcommand missing or empty!');
end

if isempty(rpixstatus) || (rpixstatus.refcount == 0)
    % Driver closed at startup:
    rpixstatus.refcount = 0;
    rpixstatus.pendingForFlip = -1;
    rpixstatus.numSamples = 1000;
    rpixstatus.bufferAddress = 12e6;
    rpixstatus.nrButtons = 5;
end

if strcmpi(cmd, 'Open')
    % Open the box: Performs basic init and setup.
    if rpixstatus.refcount > 0
        error('Open: Tried to open ResponsePixx, but ResponsePixx already open!');
    end
    
    % Open device connection if not already open:
    PsychDataPixx('Open');
    
    if PsychDataPixx('IsBusy')
        error('Open: Tried to open ResponsePixx, but Datapixx is busy! Screen flip pending?');
    end
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        if ~isscalar(varargin{1}) || ~isnumeric(varargin{1}) || varargin{1} < 1
            error('Open: numSamples must be a positive integer value.');
        end
        rpixstatus.numSamples = ceil(varargin{1});
    else
        rpixstatus.numSamples = 1000;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        if ~isscalar(varargin{2}) || ~isnumeric(varargin{2}) || varargin{2} < 0
            error('Open: bufferBaseAddress must be a positive integer value.');
        end
        rpixstatus.bufferAddress = ceil(varargin{2});
    else
        rpixstatus.bufferAddress = 12e6;
    end

    if length(varargin) >= 3 && ~isempty(varargin{3})
        if ~isscalar(varargin{3}) || ~isnumeric(varargin{3})
            error('Open: nrButtons must be a integer value.');
        end
        rpixstatus.nrButtons = ceil(varargin{3});
    else
        rpixstatus.nrButtons = 5;
    end
    
    if rpixstatus.nrButtons < 1 || rpixstatus.nrButtons > 16
        error('Open: Invalid number of buttons specified! Not in valid range 1 to 16.');
    end
    
    % Configure digital input system for monitoring the button box
    Datapixx('DisableDoutDinLoopback');                     % Disable loopback:
    Datapixx('SetDinDataDirection', hex2dec('1F0000'));     % Drive up to 5 button lights.
    Datapixx('SetDinDataOut', 0);                           % Switch all lights off.
    Datapixx('SetDinDataOutStrength', 1);                   % Set "on" lights to max intensity.
    Datapixx('EnableDinDebounce');                          % Debounce button presses.

    % Log button presses to buffer.
    Datapixx('SetDinLog', rpixstatus.bufferAddress, rpixstatus.numSamples);
    Datapixx('StopDinLog');
    
    % Execute:
    Datapixx('RegWrRd');
    
    % Mark us as open:
    rpixstatus.refcount = rpixstatus.refcount + 1;
    
    return;    
end

if strcmpi(cmd, 'Close')
    if rpixstatus.refcount <= 0
        % Noop:
        rpixstatus.refcount = 0;
        return;
    end
    
    if PsychDataPixx('IsBusy')
        error('Close: Tried to close ResponsePixx, but Datapixx is busy! Screen flip pending?');
    end
    
    % Switch all button lights off:
    Datapixx('StopDinLog');
    Datapixx('SetDinDataOut', 0);
    Datapixx('SetDinDataDirection', 0);
    Datapixx('DisableDinDebounce');
    
    % Execute:
    Datapixx('RegWrRd');
    
    % Close our connection:
    PsychDataPixx('Close');
    
    % Closed:
    rpixstatus.refcount = 0;
    
    return;
end

% All following commands need the ResponsePixx to be open:
if rpixstatus.refcount <= 0
    error('Tried to use ResponsePixx, but ResponsePixx not open! Call ResponsePixx(''Open'') first.');
end

if strcmpi(cmd, 'StartNow') || strcmpi(cmd, 'StopNow')
    if strcmpi(cmd, 'StartNow')
        startIt = 1;
    else
        startIt = 0;
    end
    
    if length(varargin) >= 1 && ~isempty(varargin{1})
        if ~isscalar(varargin{1}) || (~isnumeric(varargin{1}) && ~islogical(varargin{1}))
            error('%s: clearLog must be a integer or logical value.', cmd);
        end        
        clearLog = varargin{1};
    else
        % Don't clear log by default:
        clearLog =0;
    end

    if length(varargin) >= 2 && ~isempty(varargin{2})
        if ~isvector(varargin{2}) || ~isnumeric(varargin{2})
            error('%s: buttonLightState must be a integer vector.', cmd);
        end        
        buttonLightState = varargin{2};
    else
        buttonLightState = [];
    end

    if length(varargin) >= 3 && ~isempty(varargin{3})
        if ~isscalar(varargin{3}) || ~isnumeric(varargin{3}) || varargin{3} < 0 || varargin{3} > 1
            error('%s: buttonLightIntensity must be a value between 0.0 and 1.0.', cmd);
        end        
        buttonLightIntensity = varargin{3};
    else
        buttonLightIntensity = [];
    end

    if PsychDataPixx('IsBusy')
        error('%s: Datapixx is busy! Screen flip pending?', cmd);
    end
    
    % Build command string for setup and start/stop of logging:
    cmdstring = buildRPixxCommand(rpixstatus, startIt, clearLog, buttonLightState, buttonLightIntensity);
    
    % Execute command string:
    eval(cmdstring);
    
    % Commit:
    Datapixx('RegWrRd');
    
    if nargout > 0
        % Return boxtime of start/stop of logging:
        tStartBox = Datapixx('GetTime');
        varargout{1} = tStartBox;
    end
    
    if nargout > 1
        % Return GetSecs mapped time as well:
        tStartGetSecs = PsychDataPixx('FastBoxsecsToGetsecs', tStartBox);
        varargout{2} = tStartGetSecs;
    end
    
    if startIt
        % Mark logging as started:
        rpixstatus.pendingForFlip = 0;
    else
        % Mark logging as stopped:
        rpixstatus.pendingForFlip = -1;
    end
    
    return;
end

if strcmpi(cmd, 'StartAtFlip') || strcmpi(cmd, 'StopAtFlip')
    if strcmpi(cmd, 'StartAtFlip')
        startIt = 1;
    else
        startIt = 0;
    end
    
    if length(varargin) < 1 || isempty(varargin{1}) || varargin{1} == 0
        targetFlip = [];
    else
        if ~isscalar(varargin{1}) || ~isnumeric(varargin{1}) || varargin{1} < 0
            error('%s: flipCount must be a positive integer value.', cmd);
        end        
        targetFlip = ceil(varargin{1});
    end
    
    if targetFlip <= PsychDataPixx('FlipCount')
        % targetFlip set to a count that is already over. Pointless! Set it
        % to "next flip" instead:
        targetFlip = [];
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        if ~isscalar(varargin{2}) || (~isnumeric(varargin{2}) && ~islogical(varargin{2}))
            error('%s: clearLog must be a integer or logical value.', cmd);
        end        
        clearLog = varargin{2};
    else
        % Don't clear log by default:
        clearLog =0;
    end

    if length(varargin) >= 3 && ~isempty(varargin{3})
        if ~isvector(varargin{3}) || ~isnumeric(varargin{3})
            error('%s: buttonLightState must be a integer vector.', cmd);
        end        
        buttonLightState = varargin{3};
    else
        buttonLightState = [];
    end

    if length(varargin) >= 4 && ~isempty(varargin{4})
        if ~isscalar(varargin{4}) || ~isnumeric(varargin{4}) || varargin{4} < 0 || varargin{4} > 1
            error('%s: buttonLightIntensity must be a value between 0.0 and 1.0.', cmd);
        end
        buttonLightIntensity = varargin{4};
    else
        buttonLightIntensity = [];
    end

    if PsychDataPixx('IsBusy')
        error('%s: Datapixx is busy! Screen flip pending?', cmd);
    end
    
    % Build command string for setup and start/stop of logging:
    cmdstring = buildRPixxCommand(rpixstatus, startIt, clearLog, buttonLightState, buttonLightIntensity);
    
    % Schedule cmdstring for execution at stimulus onset for the 'targetFlip'th
    % flip command execution:
    PsychDataPixx('ExecuteAtFlipCount', targetFlip, cmdstring);
    rpixstatus.pendingForFlip = targetFlip;
    
    % Ready.
    return;
end

if strcmpi(cmd, 'GetLoggedResponses')
    if length(varargin) < 2 || isempty(varargin{2})
        blocking = 1;
    else
        if ~isscalar(varargin{2}) || (~isnumeric(varargin{2}) && ~islogical(varargin{2}))
            error('%s: blocking must be a integer or logical value.', cmd);
        end        
        blocking = varargin{2};
    end

    if length(varargin) < 3 || isempty(varargin{3})
        timeout = inf;
    else
        if ~isscalar(varargin{3}) || ~isnumeric(varargin{3}) || varargin{3} < 0
            error('%s: timeout must be a positive value.', cmd);
        end
        timeout = varargin{3};
    end
    
    if PsychDataPixx('IsBusy')
        error('%s: Datapixx is busy! Screen flip pending?', cmd);
    end
        
    % Ok, should be already running or it is at least certain that it will
    % eventually run due to trigger by psync mechanism.
    %
    % Enter a polling loop to wait until the requested minimum amount of
    % data is available:
    if length(varargin) >= 1 && ~isempty(varargin{1})
        % Wait or poll for requested amount of data:
        numFrames = ceil(varargin{1});
        if ~isscalar(numFrames) || ~isnumeric(numFrames) || numFrames < 0
            error('%s: numberResponses must be a positive value.', cmd);
        end
        
        timeout = GetSecs + timeout;
        
        while 1
            Datapixx('RegWrRd');
            status = Datapixx('GetDinStatus');
            if status.newLogFrames >= numFrames
                % Requested amount of data available: Exit loop.
                break;
            else
                % Insufficient amount. If this is a polling request, we
                % simply return no result:
                if ~blocking || (GetSecs > timeout)
                    varargout{1} = [];
                    varargout{2} = [];
                    varargout{3} = [];
                    return;
                end
            end
            
            % Sleep a msec, then retry:
            WaitSecs('YieldSecs', 0.001);
        end
    else
        % Just get current status - current amount of audiodata:
        Datapixx('RegWrRd');
        status = Datapixx('GetDinStatus');
        numFrames = -1;
    end
    
    % Either requested amount of data available or no specific amount
    % requested.

    % Fetch all captured frames if no specific amount requested:
    if numFrames == -1
        numFrames = status.newLogFrames;
    end
    
    % Retrieve logged data:
    [data, timestamps, underflow] = Datapixx('ReadDinLog', numFrames);
    
    % Decode into button state vector:
    buttons = zeros(length(data), rpixstatus.nrButtons);
    for i=1:length(data)
        for j=0:rpixstatus.nrButtons-1
            buttons(i, j+1) = ~bitand(data(i), 2^j);
        end
    end
    
    % Return all data:
    varargout{1} = buttons;
    varargout{2} = timestamps;
    varargout{3} = underflow;
    
    if underflow > 0
        fprintf('ResponsePixx: WARNING! Undeflow of event buffer detected during call to %s. Button events may be lost or corrupted!\n', cmd);
    end
    
    return;
end

if strcmpi(cmd, 'GetButtons')
    if PsychDataPixx('IsBusy')
        error('%s: Datapixx is busy! Screen flip pending?', cmd);
    end

    % Query!
    Datapixx('RegWrRd');
    
    % Retrieve current button state:
    data = Datapixx('GetDinValues');

    % Decode into button state vector:
    buttons = zeros(1, rpixstatus.nrButtons);
    for j=0:rpixstatus.nrButtons-1
        buttons(1, j+1) = ~bitand(data, 2^j);
    end
    varargout{1} = buttons;

    if nargout > 1
        % Return boxtime of query:
        tStartBox = Datapixx('GetTime');
        varargout{2} = tStartBox;
    end
    
    if nargout > 2
        % Return GetSecs mapped time as well:
        tStartGetSecs = PsychDataPixx('FastBoxsecsToGetsecs', tStartBox);
        varargout{3} = tStartGetSecs;
    end
    
    return;
end

error('ResponsePixx: Unknown subcommand provided!');

% End of main function:
end

function cmdStr = buildRPixxCommand(rpixstatus, startLogging, clearLog, buttonLightState, buttonLightIntensity)
    cmdStr = '';
    
    if ~isempty(clearLog) && (clearLog == 1)
        % Add a setup command to (re-)alloc and clear buffer:
        cmdStr = sprintf('Datapixx(''SetDinLog'', %i, %i); ', rpixstatus.bufferAddress, rpixstatus.numSamples);
    end
    
    if ~isempty(buttonLightState)
        % Setup output state for button lights:
        
        % Passing a shorter state vector than the number of buttons is
        % allowed, we will set all unspecified button lights to off. Also
        % passing a longer vector is allowed - we clamp to actual number of
        % buttons:
        nrButtons = min(length(buttonLightState), rpixstatus.nrButtons);
        
        v = 0;
        for j=1:nrButtons
            v = v + (buttonLightState(j) * 2^(15+j));
        end
        cmdStr = [ cmdStr sprintf('Datapixx(''SetDinDataOut'', %i); ', v) ];
    end
    
    if ~isempty(buttonLightIntensity)
        % Setup output intensity for button lights:
        cmdStr = [ cmdStr sprintf('Datapixx(''SetDinDataOutStrength'', %f); ', buttonLightIntensity) ];
    end

    if startLogging
        cmdStr = [ cmdStr 'Datapixx(''StartDinLog''); '];
    else
        cmdStr = [ cmdStr 'Datapixx(''StopDinLog''); '];
    end
    
    % cmdStr ready for eval()uation.
    return;
end
