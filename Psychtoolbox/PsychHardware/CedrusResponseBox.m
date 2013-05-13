function varargout = CedrusResponseBox(cmd, varargin)
% CedrusResponseBox - Interface to Cedrus Response Boxes.
%
% This function provides an interface to response button boxes from Cedrus,
% specifically model RB 530,...,830 and compatible models supporting the
% XID protocol (see http://www.cedrus.com).
%
% These response boxes connect via a serial line link interface, or a USB
% interface which emulates a serial link interface. They support the XID
% protocol for communication. See http://www.cedrus.com/xid for details.
%
% This function allows to establish a connection to the box, control a few
% of its parameters and most importantly query its button state and
% associated button press timestamps.
%
% It supports multiple subcommands, which accept and return different
% arguments, as listed below.
%
% Limitations:
% ------------
%
% Functionality is currently limited mostly to button queries (and RJ-45
% connector state queries), including timestamps, as well as control of
% built-in timers of the box. We also support basic configuration of TTL
% ports, but not yet all settings of the box like e.g., button debounce
% time. Adding such calls is straightforward and simple.
%
% We found communication with the Cedrus boxes to be unreliable quite
% often. It is an open question if this is a flaw in the design of the
% Cedrus devices and their firmware or protocols, or if the programming
% documentation for them is incomplete and therefore our implementation of
% the driver. However, the problems were reproduced under different
% operating systems, serial port drivers, toolboxes by different
% implementations written by different people, so it doesn't seem to be a
% simple glitch in one implementation. In general, the boxes work, but
% don't be surprised if you need to restart your script multiple times
% before you can establish communication, or if the more advanced
% fucntions, e.g., for configuration of the TTL RJ-45 connector, work
% unreliably for no apparent reason. Cedrus has been contacted, but so far
% no resolution or response from them.
%
% In short: If you are looking for a reliable response box that is painfree
% to use, don't buy Cedrus devices!
%
%
% Subfunctions and their meaning:
% -------------------------------
%
% Functions for device init and shutdown: Call once at beginning/end of
% your script. These are slow!
%
% handle = CedrusResponseBox('Open', port [, lowbaudrate]);
% - Open a compatible response box which is connected to the given named
% serial 'port'. 'port'names differ accross operating systems. A typical
% port name for Windows would be 'COM2', whereas a typical port name on OS/X
% or Linux would be the name of a serial port device file, e.g.,
% '/dev/cu.usbserial-FTDI125ZX9' on OS/X, or '/dev/ttyS0' on Linux.
%
% All names on OS/X are like '/dev/cu.XXXXX', where the XXXXX part depends
% on your serial port device, typically '/dev/cu.usbserial-XXXXX' for
% serial over USB devices with product name XXXXX.
%
% On Linux, all names are of pattern '/dev/ttySxx' for standard serial
% ports, e.g., '/dev/ttyS0' for the first serial port in the system, and of
% type '/dev/ttyUSBxx' for serial over USB devices, e.g., '/dev/ttyUSB0'
% for the first serial line emulated over the USB protocol.
%
%
% After the connection is established and some testing and initialization is,
% done, the function returns a device 'handle', a unique identifier to use
% for all other subfunctions.
%
% By default the commlink is opened at a baud transmission rate of 115200
% Baud (All DIP switches on the box need to be in 'down' position!). If you
% specify the optional flag 'lowbaudrate' as 1, then the speed will be
% lowered to 56 kBaud at device open time -- in case your system works
% unreliably at the higher rate.
%
% By default, the script uses Psychtoolbox's own IOPort() serial link
% driver for communication (ptb_cedrus_drivertype = 2). If you want to use
% a different driver for testing, change the 'ptb_cedrus_drivertype'
% parameter inside the code with the id of a supported driver (Matlab
% serial() on Windows and Linux, SerialComm on OS/X). This option may go
% away in the future and is for debugging only!
%
%
% CedrusResponseBox('Close', handle);
% - Close connection to response box. The 'handle' becomes invalid after
% that command.
%
%
% CedrusResponseBox('CloseAll');
% - Close all connections to all response boxes. This is a convenience
% function for quick shutdown.
%
%
% dev = CedrusResponseBox('GetDeviceInfo', handle);
% - Return queried information about the device in a struct 'dev'. 'dev'
% contains (amongst other) the following fields:
%
% General information:
%  dev.Name = Device name string.
%  dev.VersionMajor and dev.VersionMinor = Major and Minor firmware revision.
%  dev.productId = Type of device, e.g., 'Lumina', 'VoiceKey' or 'RB response pad'.
%  dev.modelId   = Submodel of the device if the device is a RB response pad,
%                  e.g., 'RB-530', 'RB-730', 'RB-830' or 'RB-834'.
%
%  dev.port      = Portname of serial port, as passed to the open function.
%
% Diagnostic information for timing: Values of -1 or 0 usually mean "info
% not available".
%
%  dev.roundtriptime   = Median of estimated roundtrip latency for
%                        communication with the box - in seconds.
%
%  dev.roundtripstddev = Standard deviation from mean of roundtrip latency
%  measurements in seconds. Large numbers mean that your operating system
%  has bad scheduling and that reported event timestamps may be uncertain by
%  that amount.
%
%  dev.rttresetdelay   = Duration (in seconds) of a reaction time timer reset sequence
%  Values of more than 3 msecs indicate some problems with the box itself or
%  the communication link -- Measured event times or reaction times may not
%  be trustworthy!
%
%
% Functions for use within script. These are as fast as possible:
%
% CedrusResponseBox('ClearQueues', handle);
% - Clear all queues, discard all pending data.
%
% [status = ] CedrusResponseBox('FlushEvents', handle);
% - Empty/clear/flush the queue of pending events. Use this to get rid of
% any stale button press or release events before start of response
% collection in a trial. E.g., Assume you wait for a subjects keypress and
% finally receive that keypress via 'GetButtons' or 'WaitButtons'. You
% collected your response, the trial is done, but when the subject releases
% the button again, that will generate another event - a release event, in
% which you're not interested. Maybe the subject will accidentally hit the
% button as well. --> Good to clean the queue before a new trial.
%
% This function has a second use as well. It has an optional output
% argument, 'status', which will return the current status of all buttons
% (i.e. whether they are currently being pressed or not).
% Status is a 3 row by 8 column matrix: Row 1 describes the status of the
% up to eight pushbuttons of the box. Row 2 describes the status of the TTL
% lines of the RJ-45 accessory connector. Row 3 describes the status of the
% VoiceKey if any. Columns 1 to 8 of each row correspond to buttons 1-8,
% TTL lines 1-8 or inputs 1-8 of the VoiceKey.
%
% The mapping for the CB-530 for row 1 of 'status' status(1,:) is as follows:
%
% [top ??? left middle right bottom] -- the 2nd entry has no associated
% button, but it may be the scanner trigger input. The mapping on other boxes
% may be different.
%
% This is useful if you just want to know whether the subject is currently
% pressing any buttons before you proceed, but are not fussed about timing.
%
% E.g. I often find myself doing the following:
%   buttons = 1;
%   while any(buttons(1,:))
%     buttons = CedrusResponseBox('FlushEvents', mybox);
%   end
%
% ...to wait for the subject to release any buttons which might currently be down.
%
% evt = CedrusResponseBox('GetButtons', handle);
% - Return next queued button-press or button-release event from the box.
% Each time a button on the box is pressed or released, and each time the
% state of the accessory connector changes, an "event" data packet is sent
% from the box to the computer. The packet is timestamped with the time of
% the triggering event, as measured by the boxes reaction time timer.
%
% This function checks if such an event is available and returns its
% description in a 'evt' struct, if so. If no event is pending, it returns an
% empty 'evt', ie. isempty(evt) is true.
%
% 'evt' for a real fetched event is a struct with the following fields:
%
% evt.raw     = "raw" byte that describes the event. Only for debugging.
%
% evt.port    = Number of the device port on which the event occured. Push
%               buttons and scanner triggers are on port 0, the RJ-45 TTL
%               connector is on port 1, port 2 is the voice-key (if any).
%
% evt.action  = Action that triggered the event:
%               1 = Button press, 0 = Button release for pushbuttons.
%               1 = TTL line high, 0 = TTL line low for RJ-45 I/O lines.
%               1 = Voice onse, 0 = Voice offset/silence for Voicekey.
%
% evt.button  = Number of the button that was pressed or released (1 to 8)
%               or the TTL line that was going high/low. Numbers vary by
%               response box.
%
% evt.buttonID= Descriptive name string for pressed button, e.g., 'top' or
%               'left'. Please note that this mapping is only meaningful
%               for the RB-530 response box.
%
% evt.rawtime = Time of the event in secs since last reset of the reaction
%               time timer, measured in msecs resolution. This value is
%               always valid, but not directly comparable to any other
%               timestamps or time measurements within Psychtoolbox.
%
%
% evt = CedrusResponseBox('WaitButtons', handle);
% - Queries and returns the same info as 'GetButtons', but waits for
% events. If there isn't any event available, will wait until one becomes
% available.
%
% evt = CedrusResponseBox('WaitButtonPress', handle);
% - Like WaitButtons, but will wait until the subject /presses/ a key -- the
% signal that a key has been released is not acceptable -- Button release
% events are simply discarded.
%
%
% evt = CedrusResponseBox('GetBaseTimer', handle [, nSamples=1]);
% - Query current time of base timer of the box. Returned values are in
% seconds, resolution is milliseconds. evt.basetimer is the timers time,
% maybe corrected for serial link receive latency. evt.ptbreceivetime is a
% timestamp taken via PTB's GetSecs() at time of receive of the data.
% evt.ptbtime is the basetimers time mapped into PTB GetSecs time if such a
% mapping is possible, otherwise this field doesn't exist:
% evt.ptbreceivetime and evt.ptbtime shouldn't be significantly different
% if everything is good. Large differences indicate some timing problems
% with the connection to the box, or a timer problem - either with your
% computers timer or the hardware timer of the tox, or significant
% clock-drift between the computers timer and the boxes timer. In any case,
% reaction timer measurements and such will be problematic.
%
% Note that this automatically discards all pending events in the queue before
% performing the timer query!
%
% The optional argument 'nSamples' allows to specify if multiple samples of
% PTB timer vs. the response boxes timer should be measured. If 'nSamples'
% is set to a value greater than one, a cell array with nSamples elements
% will be returned, each corresponding to one measurement. This allows,
% e.g., to check if PTBs timer and the boxes timer drift against each
% other.
%
%
% resetTime = CedrusResponseBox('ResetRTTimer', handle);
% - Reset reaction time timer of box to zero. This should not be neccessary
% if you use the evt.ptbtime timestamps for time measurements or reaction
% time measurements. If you however use uncalibrated mode and the
% evt.rawtime values directly, this function may be useful to establish a
% zero baseline for reaction time measurements. However, as the communication
% delay for sending the reset command can't be reliably measured, using
% such a software triggered timer reset may not be the most reliable way of
% resetting the timer. The function returns 'resetTime' PTB's best guess of
% when the reset was carried out -- essentially a GetSecs() timestamp of
% when the reset command was sent.
%
% Note that this automatically discards all pending
% events in the queue before performing the query!
%
%
% slope = CedrusResponseBox('GetBoxTimerSlope', handle);
% - Compute slope (drift) between computer clock and device clock. 'slope'
% tells how many seconds of time "elapse" on the computer in GetSecs time
% for each "elapsed" second of box time. At device open time, the driver
% takes a timestamp from the device basetimer. This function also takes a
% timestamp and then computes the ratio of differences. The longer you'll
% wait after CedrusResponseBox('Open') before calling this function, the
% more accurate the clock-drift estimate will be.
%
%
% roundtrip = CedrusResponseBox('RoundTripTest', handle);
% - Initiate 100 trials of the roundtrip test of the box. Data is echoed
% forth and back 100 times between PTB and the box, and the latency is
% measured (in seconds, with msecs resolution). The vector of all samples
% is returned in 'roundtrip' for evaluation and debugging. The measured
% latency is also used for delay correction for the 'GetBaseTimer'
% subfunction. However, a roundtrip test is performed automatically when
% opening the response box connection, so this is rarely needed.
%
% Note that this automatically discards all pending
% events in the queue before performing the query!
%
%
% [currentMode] = CedrusResponseBox('SetConnectorMode', handle [, mode]);
% - Set or get mode of operation of external accessory connector: 'mode' can be
% any of the following text strings:
%
% 'GeneralPurpose': Input/Output assignment of pins can be freely
% programmed via the 'DefineInputLinesAndLevels' subcommand (see below),
% and the output lines only change if the 'SetOutputLineLevels' command
% (see below) is used. The connector doesn't change state by itself.
%
% 'ReflectiveContinuous': Line levels reflect button state: Line is active
% if button is pressed and goes inactive when the button is released again.
%
% 'ReflectiveSinglePulse': A single pulse is sent to an output line if a
% button is pressed on the box. Nothing is sent on release.
%
% 'ReflectiveDoublePulse': A single pulse is sent to an output line if a
% button is pressed on the box. Another pulse is sent on button release.
%
% If 'mode' is left out, the function queries and returns the current mode
% as return argument 'currentMode'. If mode is given, nothing is returned.
%
%
% CedrusResponseBox('SetOutputLineLevels', handle, outlevels);
% - Set accessory connector output lines to state specified in 'outlevels'.
% outlevels is an 8 element vector of zeros and ones. Each element
% corresponds to an output pin, and its values sets the output level of
% that pin. Example: outlevel = [1,1,1,1,0,0,0,0] would set the 4 lines
% with the lowest numbers (lines 0,1,2,3) to '1' aka active and the 4 lines
% with the highest numbers (lines 4,5,6,7) to '0' aka inactive.
% This corresponds to XiD command 'ah'.
%
% The command is only effective if connector is set to 'GeneralPurpose'.
%
%
% CedrusResponseBox('DefineInputLinesAndLevels', handle, inputlines, logiclevel, debouncetime);
% - Define which lines on the connector are inputs: 'inputlines' is a
% vector with the line numbers of the input lines. All other lines are
% designated as output lines, e.g., inputlines = [0, 2, 4] would set lines
% 0, 2 and 4 as inputs, remaining lines 1,3,5,6,7 as outputs. 'logiclevel'
% tells if the default TTL level of the input lines is low (logiclevel=1)
% or high (logiclevel=0). Example: logiclevel = 1 means that the lines are
% pulled low by default, so they will detect an active high state -- if
% their level is raised to TTL high state. The argument 'debouncetime' must
% be the debounce time for the input lines in milliseconds. After an event
% on a input line, the box will ignore all further events on than input
% line for 'debouncetime' milliseconds.
%
% This corresponds to XiD commands 'a4', 'a50' and 'a51', as well as 'a6'.
%
% The command is only effective if connector is set to 'GeneralPurpose'.
%
%
% inputLines = CedrusResponseBox('ReadInputLines', handle);
% - Read current state of the connectors input lines: Returns an 8 element
% vector where each element corresponds to one input line and a 1 means
% active, 0 means inactive. This corresponds to XiD command 'ar'.
%
% Note that this automatically discards all pending
% events in the queue before performing the query!
%
% The command is only effective if connector is set to 'GeneralPurpose'.
%
%

% Technical notes:
% USB VendorID of Cedrus:   0x0403
% USB ProductID:            0xf228
%
% Command for manual insertion of serial-over-USB module on Linux, if
% module doesn't recognize Cedrus device id's. Also edit rules files of
% usbdev, so ftdi_sio module gets auto-loaded on Cedrus insertion.
%
% sudo modprobe ftdi_sio product=0xf228

% Disabled help text snippets:
% [,doCalibrate=0]
% If you don't specify the optional 'doCalibrate' flag, or
% leave it at its default setting of 1, a couple of lengthy (multiple
% seconds) timing calibrations and tests are performed. These allow to
% assess the delays in communication between box and Matlab. They will also
% allow to return all times of events (as detected by the box) in PTB's
% standard GetSecs() time reference system -- Timestamps of button press
% events and TTL input events can be directly compared with timestamps
% delivered by other PTB functions like GetSecs, KbCheck, KbWait,
% Screen('Flip') etc.
%
% If you set the 'doCalibrate' flag to zero, all timing calibrations will
% be skipped: Startup time is drastically reduced. However there isn't any
% simple and straightforward way of comparing timestamps or timer readings
% delivered by the box with other timestamps of PTB functions. This only
% makes sense if you use some external triggering mechanism to reset the
% built-in reaction time timer via some external TTL input trigger signals
% and want to use raw timer measurements.

% evt.ptbtime = Time of the event in secs, measured in PTBs "GetSecs"
%               timebase. This is easier to correlate with other
%               timestamps, e.g., Screen('Flip') timestamps, but its
%               reliability hasn't been tested yet for the current
%               software release. When opening a connection to a response
%               box, we perform timing calibrations to establish the
%               mapping of time values as measured by the hardware timers
%               of your response pad to time values in PTB's reference
%               system. If you skipped that calibrations by setting the
%               optional 'doCalibrate' flag to zero at device open time,
%               then the evt.ptbtime field will not be available and you
%               have to cope with evt.rawtime values only.
%


% History:
%
% 03/21/08 Written. Based on example code donated by Cambridge Research Systems. (MK)
% 03/28/08 Altered by Jenny Read.
% 04/03/08 Refined and added MacOS/X support via SerialComm driver. (MK)
% 04/06/08 Improved timing code for mapping of box timers --> GetSecs time. (MK)
% 04/17/08 Disable Boxtime->Ptbtime mapping for now, use old drivers. (MK)
% 04/23/08 Add additional setup and query commands for external port. (MK)
% 05/09/10 Add additional button label definitions for RB830, contributed
%          by Jochen Laubrock. (MK)

% Hard-Coded drivertype to use: Defaults to our IOPort driver.
global ptb_cedrus_drivertype;

ptb_cedrus_drivertype = 2;

% Cell array of device structs. Globally available for main function and
% all subfunctions in this file, persistent across invocation:
global ptb_cedrus_devices;

% Subcommand dispatch:
if nargin < 1 || ~ischar(cmd)
    error('Must at least specify subcommand as textstring!')
end

% Following if-end blocks are roughly sorted by frequency of use. The most
% frequently used calls come first to achieve minimum dispatch
% latency in the trial loop.

if strcmpi(cmd, 'FlushEvents')
    % Flush all pending events/data:

    if nargin < 2
        error('You must provide the device "handle" for the box to flush!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Perform flush:
    CedrusStatus = FlushEvents(handle);

    if nargout>0
        varargout{1} = CedrusStatus;
    end
    return
end

if strcmpi(cmd, 'ClearQueues')
    % Clear all pending events/data:

    if nargin < 2
        error('You must provide the device "handle" for the box to clear!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Perform clear:
    ClearQueues(handle);

    return
end

% Wait until a key-pressed signal is detected:
if strcmpi(cmd, 'WaitButtonPress')
    if nargin < 2
        error('You must provide the device "handle" for the box to wait for!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Wait until the box reports that a key has been pressed (NOT released)
    keypress = 0;
    while ~keypress
        % Need at least 6 bytes. Only wait if not available:
        if BytesAvailable(handle) < 6
            % Poll at 4 msecs intervals as long as input buffer is totally empty,
            % to allow the CPU to execute other tasks.
            while BytesAvailable(handle) == 0
                % Choose 4 msecs, as PTB would not release the cpu for wait
                % times below 3 msecs (to account for MS-Windows miserable
                % process scheduler).
                if IsWin
                    WaitSecs(0.004);
                else
                    WaitSecs(0.001);
                end
            end

            % At least 1 byte available -- soon we'll have our required minimimum 6
            % bytes :-) -- Spin-Wait for the remaining few microseconds:
            while BytesAvailable(handle) < 6; end
        end

        % At least 6 bytes for one event available: Try to read them from box:
        response = ReadDev(handle, 6);
        
        % Timestamp receive completion in PTB's timeframe. Allows to get a
        % feeling on how much time elapses between keypress and data receive:
        ptbfetchtime = GetSecs;

        % Unpack this binary data into a more readable form:
        evt = ExtractKeyPressData(handle,response);        
        evt.ptbfetchtime = ptbfetchtime;
        
        keypress = evt.action;
        % This is 0 if the key was released, 1 if it was pressed down,
        % which is what we are waiting for.
    end
    % Assign evt as output argument:
    varargout{1} = evt;

    return % JCAR added
end

% Wait for at least one button event available:
if strcmpi(cmd, 'WaitButtons')
    if nargin < 2
        error('You must provide the device "handle" for the box to wait for!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Need at least 6 bytes. Only wait if not available:
    if BytesAvailable(handle) < 6
        % Poll at 4 msecs intervals as long as input buffer is totally empty,
        % to allow the CPU to execute other tasks.
        while BytesAvailable(handle) == 0
            % Choose 4 msecs, as PTB would not release the cpu for wait
            % times below 3 msecs (to account for MS-Windows miserable
            % process scheduler).
            if IsWin
                WaitSecs(0.004);
            else
                WaitSecs(0.001);
            end
        end

        % At least 1 byte available -- soon we'll have our required minimimum 6
        % bytes :-) -- Spin-Wait for the remaining few microseconds:
        while BytesAvailable(handle) < 6; end
    end

    % At least 6 bytes for one event available: Try to read them from box:
    response = ReadDev(handle, 6);

    % Timestamp receive completion in PTB's timeframe. Allows to get a
    % feeling on how much time elapses between keypress and data receive:
    ptbfetchtime = GetSecs;

    % Unpack this binary data into a more readable form:
    evt = ExtractKeyPressData(handle,response);
    evt.ptbfetchtime = ptbfetchtime;

    % Assign evt as output argument:
    varargout{1} = evt;

    return % JCAR added
end

% Polling Button state query: Returns immediately if no events available.
if strcmpi(cmd, 'GetButtons')
    % Button state change event query:

    if nargin < 2
        error('You must provide the device "handle" for the box to query!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % XID devices send six bytes of information.
    % Check if at least 6 bytes for a full event report are available:
    if BytesAvailable(handle) < 6
        % Did not receive at least 6 bytes - No new event happened. Return an
        % empty evt result.
        evt = [];
        varargout{1} = evt;
        return;
    end

    % At least 6 bytes for one event available: Try to read them from box:
    response = ReadDev(handle, 6);

    % Timestamp receive completion in PTB's timeframe. Allows to get a
    % feeling on how much time elapses between keypress and data receive:
    ptbfetchtime = GetSecs;
    
    % Unpack this binary data into a more readable form:
    evt = ExtractKeyPressData(handle,response);
    evt.ptbfetchtime = ptbfetchtime;

    % Assign evt as output argument:
    varargout{1} = evt;

    return;
end

if strcmpi(cmd, 'RoundTripTest')
    % Initiate roundtrip-test procedure: Will receive data from device,
    % echo it back, then receive a roundtrip timestamp:

    if nargin < 2
        error('You must provide the device "handle" for the box to query!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});
        
    % Start roundtrip test and return results: They are also stored in the
    % device struct of 'handle':
    varargout{1} = RoundTripTestDev(handle);
    return;
end

if strcmpi(cmd, 'ResetRTTimer')
    % RT timer reset request:
    if nargin < 2
        error('You must provide the device "handle" for the box!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Reset reaction time timer of device and assign estimated time of reset
    % as basetime for all timing calculations:
    varargout{1} = ResetRTT(handle);

    return;
end

if strcmpi(cmd, 'Test')
    % Flush all pending events/data:

    if nargin < 2
        error('You must provide the device "handle" for the box!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    TestThis(handle);
    
    return
end

if strcmpi(cmd, 'GetBaseTimer')
    % Base Timer query:

    if nargin < 2
        error('You must provide the device "handle" for the box to query!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    if nargin >=3
        nQueries = varargin{2};
    else
        nQueries = 1;
    end
    
    % Preallocate output cell array:
    evts = cell(nQueries, 1);
    
    % Flush input buffer:
    FlushEvents(handle);

    for i=1:nQueries
        % Send basetimer query code:
        evt.roundtriptime = GetSecs;
        WriteDev(handle, 'e3');

        % Spin-Wait for first byte:
        while BytesAvailable(handle) < 1; end;

        % Timestamp receive completion of first byte. This is closest to the
        % real time when the transmitted timer values was actually generated on
        % the device:
        evt.ptbreceivetime = GetSecs;

        % Receive packet, then parse into raw timer value (in seconds):
        evt.basetimer = receiveAndParseTimePacket(handle);

        % Store roundtrip-time of query:
        evt.roundtriptime = evt.ptbreceivetime - evt.roundtriptime;
        
        % Correct reported time value of basetimer by half roundtrip delay
        % of serial link: We assume that transmission took half the total
        % measured roundtrip time, so we need to add that delay to the
        % basetimer value to get an estimate of the "real" basetimer time
        % at time of response packet receive "ptbtime":
        evt.basetimer = evt.basetimer + ptb_cedrus_devices{handle}.roundtriptime/2;

        % Assign mapped PTB GetSecs time if mapping possible:
        if ptb_cedrus_devices{handle}.baseToPtbSlope ~= 0
            % Simple linear equation mapping:
            evt.ptbtime = ptb_cedrus_devices{handle}.baseToPtbOffset + ptb_cedrus_devices{handle}.baseToPtbSlope * evt.basetimer;
        end
        
        % Assign i'th measurement event:
        evts{i} = evt;
    end
    
    % Assign evts as output argument:
    if nQueries > 1
        varargout{1} = evts;
    else
        varargout{1} = evt;
    end
    
    return;
end

if strcmpi(cmd, 'SetConnectorMode')
    % Change mode of external accessory connector:
    if nargin < 2
        error('You must provide the device "handle" for the box!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    if nargin < 3
        % Query instead of set:

        % Send query code:
        WriteDev(handle, '_a1');
        
        % Retrieve response:
        cc = ReadDev(handle, 4);
        WaitSecs(0.25);

        switch char(cc)
            case {'_a10'}
                rc = 'generalpurpose';
            case {'_a11'}
                rc = 'reflectivecontinuous';
            case {'_a12'}
                rc = 'reflectivesinglepulse';
            case {'_a13'}
                rc = 'reflectivedoublepulse';
            otherwise
                rc = cc;
                warning('SetConnectorMode received unknown old mode response!');
        end
        
        varargout{1} = rc;
        return;
    end

    switch lower(char(varargin{2}))
        case {'generalpurpose'}
            cc = 'a10';
        case {'reflectivecontinuous'}
            cc = 'a11';
        case {'reflectivesinglepulse'}
            cc = 'a12';
        case {'reflectivedoublepulse'}
            cc = 'a13';
        otherwise
            error('Unknown connector mode specified to SetConnectorMode.');
    end

    % Send command code:
    WriteDev(handle, cc);
    WaitSecs(0.25);
    
    return;
end

if strcmpi(cmd, 'DefineInputLinesAndLevels')
    % Change I/O assignment and default logic level of pins:
    if nargin < 2
        error('DefineInputLinesAndLevels: You must provide the device "handle" for the box!');
    end

    if nargin < 3
        error('DefineInputLinesAndLevels: You must provide the new list of integer input pin numbers!');
    end

    if nargin < 4
        error('DefineInputLinesAndLevels: You must provide the new logic detection level for the inputs: 1 for "Detect Low->High transition", 0 for "Detect High->Low" !');
    end

    if nargin < 5
        error('DefineInputLinesAndLevels: You must provide the new debounce time for TTL inputs in milliseconds !');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    ClearQueues(handle);
    
    inpins = varargin{2};
    if ~isnumeric(inpins)
        error('DefineInputLinesAndLevels: You must provide the new list of integer input pin numbers!');
    end

    if ~isempty(inpins)
        if min(inpins) < 0 || max(inpins) > 5
            error('DefineInputLinesAndLevels: Only input pin numbers between 0 and 5 are valid!');
        end
    end
    
    ipin = 0;
    for i=1:length(inpins)
        ipin = ipin + 2^(inpins(i));
    end
    
    % Send command code and mask:
    WriteDev(handle, ['a4' char(ipin)]);
    
    % Wait a bit:
    WaitSecs(1);    
    ClearQueues(handle);
    WaitSecs(1);    
    
    % Retrieve new mask:
    WriteDev(handle, '_a4');
    WaitSecs(1);
    resp = ReadDev(handle, 4);
    if length(resp)<4 || ~strcmp(char(resp(1:3)), '_a4')
        warning('DefineInputLinesAndLevels: Invalid response received from device!');
        char(resp) %#ok<NOPRT>
        varargout{1} = 0;
        return;
    else
        if resp(4)~=ipin
            warning('DefineInputLinesAndLevels: Real I/O bitmask not equal to requested one!');
            resp(4) %#ok<NOPRT>
            varargout{1} = 0;
            return;
        end
    end
    WaitSecs(1);
    
    % Send new TTL pull-level:
    if varargin{3} > 0
        % Pull lines low --> Detect lines high:
        cc = 'a50';
    else
        % Pull lines high --> Detect lines low:
        cc = 'a51';
    end
    WriteDev(handle, cc);
    WaitSecs(1);
    ClearQueues(handle);    
    WriteDev(handle, '_a5');
    WaitSecs(1);
    
    % Query pull level:
    resp = ReadDev(handle, 4);
    if length(resp) < 4
        warning('DefineInputLinesAndLevels: No response received from device!');
        varargout{1} = 0;
        return;        
    end
    
    if ~strcmp(char(resp(1:3)), '_a5')
        warning('DefineInputLinesAndLevels: Invalid response received from device!');
        char(resp) %#ok<NOPRT>
        varargout{1} = 0;
        return;        
    else
        if ~strcmp(char(resp(2:4)), cc)
            warning('DefineInputLinesAndLevels: Real TTL default not equal to requested one!');
            resp(2:4) %#ok<NOPRT>
            varargout{1} = 0;
            return;
        end
    end
    WaitSecs(1);
    
    % Send new debounce time:
    WriteDev(handle, ['a6' char(double(varargin{4}))]);
    WaitSecs(1);
    ClearQueues(handle);
    % Read it back:
    WriteDev(handle, '_a6');
    WaitSecs(1);

    % Query debounce time:
    resp = ReadDev(handle, 4);
    if length(resp) < 4
        warning('DefineInputLinesAndLevels: No response received from device!');
        varargout{1} = 0;
        return;        
    end
    
    if ~strcmp(char(resp(1:3)), '_a6')
        warning('DefineInputLinesAndLevels: Invalid response received from device!');
        char(resp);
        varargout{1} = 0;
        return;        
    else
        if resp(4) ~= double(varargin{4})
            warning('DefineInputLinesAndLevels: Real TTL debounce time not equal to requested one!');
            double(resp(4));
            varargout{1} = 0;
            return;
        end
    end
    WaitSecs(1);
    
    varargout{1} = 1;

    return;
end

if strcmpi(cmd, 'ReadInputLines')
    % Retrieve state of all input lines:
    if nargin < 2
        error('You must provide the device "handle" for the box!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Flush event queue:
    FlushEvents(handle);
    
    % Send command code:
    WriteDev(handle, 'ar');
    
    % Read one byte back:
    inplines = ReadDev(handle, 1);
    outv=zeros(1,8);
     
    for i=0:7
        if bitand(inplines, 2^i)
            outv(i+1)=1;
        end
    end
    
    varargout{1} = outv;
    
    return;
end

if strcmpi(cmd, 'SetOutputLineLevels')
    % Change signal level of output pins:
    if nargin < 2
        error('You must provide the device "handle" for the box!');
    end

    if nargin < 3
        error('You must provide the 8 element vector of output line levels!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    opins = varargin{2};
    if length(opins)~=8 || ~isnumeric(opins)
        error('You must provide an 8 element vector of output line levels!');
    end

    outval = 0;
    for i=1:8
        if opins(i)>0
            outval = outval + 2^(i-1);
        end
    end
    
    % Send command code and mask:
    WriteDev(handle, ['ah' char(outval)]);
    
    % Wait a bit:
    WaitSecs(0.1);
    
    return;
end



if strcmpi(cmd, 'GetDeviceInfo')
    % Query info about device:

    if nargin < 2
        error('You must provide the device "handle" for the box to query!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    dev = ptb_cedrus_devices{handle};
    
    varargout{1} = dev;

    return; % JCAR added
end

if strcmpi(cmd, 'Open')

    % Open a new connection to response box:

    if nargin < 2
        error('You must provide the "port" parameter for the serial port to which the box is connected!')
    end

    % Create serial object for provided port, configure connection
    % properly:
    port = varargin{1};
    % port = '/dev/cu.usbserial-FT3Z95V5'

    if nargin < 3
        % Assume user doesn't want time calibration:
        lowbaudrate = 0;
    else
        lowbaudrate = varargin{2};
    end
  
    if nargin < 4
        % Assume user doesn't want time calibration:
        doCalibrate = 0;
    else
        doCalibrate = varargin{3};
    end
    
    % Open device link at default baudrate of 115 kBaud, return 'dev' struct:
    dev = OpenDev(port, 115200);
    
    % Create new entry in our struct array:
    if isempty(ptb_cedrus_devices)
        ptb_cedrus_devices = cell(1,1);
    else
        ptb_cedrus_devices(end+1) = cell(1,1);
    end

    % Get a handle to it:
    handle = length(ptb_cedrus_devices);

    % Assign device struct to array:
    ptb_cedrus_devices{handle} = dev;
    clear dev;

    if lowbaudrate
        % Set the device protocol to XID mode
        WriteDev(handle, 'c10'); %JCAR removed cr

        % Give device time to settle:
        WaitSecs(0.5);

        % Initiate a device reset:
        WriteDev(handle, 'f7');

        % Give device time to settle:
        WaitSecs(0.5);

        % Change baudrate of device to 57600 Baud:
        WriteDev(handle, ['f1' char(3)])

        % Give device time to settle:
        WaitSecs(0.5);

        % Close connection:
        CloseDev(handle);

        % Give device time to settle:
        WaitSecs(0.5);

        % Reinit connection at new baud rate:
        % Open device link at new baudrate of 57600 Baud, return 'dev' struct:
        dev = OpenDev(port, 57600);

        % Reassign device struct to array:
        ptb_cedrus_devices{handle} = dev;
        clear dev;
    else
        % Set the device protocol to XID mode
        WriteDev(handle, 'c10'); %JCAR removed cr

        % Give device time to settle:
        WaitSecs(0.5);
        
        % Initiate a device reset:
        WriteDev(handle, 'f7');

        % Give device time to settle:
        WaitSecs(0.5);
    end
    
    % This is for keeping track of what buttons are currently up or
    % down. I assume that all buttons are up when the device is opened.
    ptb_cedrus_devices{handle}.CedrusStatus = zeros(3,8);
    
    % Put this in a try-catch loop so that if it doesn't work for any
    % reason, I can then close the link and you can try again. Otherwise,
    % the COM port is permanently busy and I have to restart Matlab.
%    try
        
        
        % % Debug information from http://www.cedrus.com/xid/properties.htm
        % %
        % % Note: 0 is ASCII value 48, 1 is ASCII value 49, and so forth.
        %
        % fprintf(s1,['_d3',char(13)]);
        % fread(s1,1)
        % fprintf(s1,['_d1',char(13)]);
        % fscanf(s1)

        % Query a few device properties:

        % Get product ID: 0 = Lumina, 1 = VoiceKey, 2 = RB response pad:

        % I have to put this in a while loop, because sometimes '_d2' fails to
        % evoke a response:
        bytes = 0;
        while bytes==0
            WriteDev(handle, '_d2');
            WaitSecs(0.25); % I also have to wait, because even when it does evoke a response,
            % there can be a long delay - tens of milliseconds. Jon Peirce confirms
            % this.
            bytes = BytesAvailable(handle);
        end

        response=ReadDev(handle, bytes);
        response=response(1);

        switch response(1)
            case 48
                ptb_cedrus_devices{handle}.productID = 'Lumina';
            case 49
                ptb_cedrus_devices{handle}.productID = 'VoiceKey';
            case 50
                ptb_cedrus_devices{handle}.productID = 'RB response pad';
            otherwise
                ptb_cedrus_devices{handle}.productID = 'Unknown';
        end

        % Get model ID: 0 = Unknown, 1 = RB-530, 2 = RB-730, 3 = RB-830, 4 = RB-834
        % I have to put this in a while loop, because sometimes '_d3' fails to
        % evoke a response:

        % Give device time to settle:
        WaitSecs(0.5);

        % Remove junk - if any:
        while BytesAvailable(handle)
            ReadDev(handle, 1);
        end
        
        bytes = 0;
        while bytes==0
            WriteDev(handle, '_d3');
            WaitSecs(0.25); % I also have to wait, because even when it does evoke a response,
            % there can be a long delay - tens of milliseconds. Jon Peirce confirms
            % same behaviour on his system.
            bytes = BytesAvailable(handle);
        end
        
        response=ReadDev(handle, bytes);
        response=response(1);

        if response==48
            ptb_cedrus_devices{handle}.modelID = 'Unknown';
            ptb_cedrus_devices{handle}.modelNo = 0;

        else if strcmp(ptb_cedrus_devices{handle}.productID,'RB response pad')
                switch response
                    case 49
                        ptb_cedrus_devices{handle}.modelID = 'RB-530';
                        ptb_cedrus_devices{handle}.modelNo = 530;
                    case 50
                        ptb_cedrus_devices{handle}.modelID = 'RB-730';
                        ptb_cedrus_devices{handle}.modelNo = 730;
                    case 51
                        ptb_cedrus_devices{handle}.modelID = 'RB-830';
                        ptb_cedrus_devices{handle}.modelNo = 830;
                    case 52
                        ptb_cedrus_devices{handle}.modelID = 'RB-834';
                        ptb_cedrus_devices{handle}.modelNo = 834;
                    otherwise
                        ptb_cedrus_devices{handle}.modelID = sprintf('Unknown id %i', response);
                        ptb_cedrus_devices{handle}.modelNo = 0;
                end
            else
                ptb_cedrus_devices{handle}.modelID = 'Unknown';
                ptb_cedrus_devices{handle}.modelNo = 0;
            end
        end

        % Firmware revision:
        bytes = 0;
        while bytes==0
            WriteDev(handle, '_d4');
            WaitSecs(0.1);
            bytes = BytesAvailable(handle);
        end
        ptb_cedrus_devices{handle}.VersionMajor = ReadDev(handle, bytes) - 48;

        bytes = 0;
        while bytes==0
            WriteDev(handle, '_d5');
            WaitSecs(0.1);
            bytes = BytesAvailable(handle);
        end
        ptb_cedrus_devices{handle}.VersionMinor = ReadDev(handle, bytes) - 48;

        % Product name string:
        bytes = 0;
        while bytes==0
            WriteDev(handle, '_d1');
            WaitSecs(0.1);
            bytes = BytesAvailable(handle);
        end
        % Weird casting procedure with replacement of char(13) by char(10),
        % so Octave can handle it:
        ptb_cedrus_devices{handle}.Name = double(ReadDev(handle, bytes));
        ptb_cedrus_devices{handle}.Name(find(ptb_cedrus_devices{handle}.Name == 13)) = 10; %#ok<FNDSB>
        ptb_cedrus_devices{handle}.Name = char(ptb_cedrus_devices{handle}.Name);
        
        % Try our best to totally drain the receive queue:
        WaitSecs(0.25);
        while 1
            bytes = BytesAvailable(handle);
            if bytes == 0
                break;
            end
            ReadDev(handle, bytes)
            WaitSecs(0.1);
        end
        
        % Reset base timer:
        WriteDev(handle, 'e1');
      
        % Calibration of PTB's timebase vs. Boxes timebase wanted?
        if doCalibrate

            % Set slope of 1 as a flag that ResetRTT should do a
            % calibrated, timestamped reset:
            ptb_cedrus_devices{handle}.baseToPtbSlope = 1;
            
            % Perform calibrated basetimer query:
            [hosttime, devicetime, minwin] = queryBaseTimer(handle);

            % And store its results:
            ptb_cedrus_devices{handle}.lastBaseTimeQuery = [hosttime, devicetime, minwin];
            
        else
            % Uncalibrated mode requested. Saves a few seconds of startup
            % time, but doesn't allow mapping of boxes time measurements
            % into GetSecs() timebase of PTB:
            
            % No link roundtrip time estimates available:
            ptb_cedrus_devices{handle}.roundtriptime = 0;
            ptb_cedrus_devices{handle}.roundtripstddev = 0;

            % No mapping of box time to PTB time available:
            ptb_cedrus_devices{handle}.baseToPtbSlope  = 0;
            ptb_cedrus_devices{handle}.baseToPtbOffset = 0;
            
            ptb_cedrus_devices{handle}.lastBaseTimeQuery = [];
        end
        
        % Reset reaction time timer of device: If calibration was
        % requested, this will also estimate the offset between RTT values
        % and basetimer values, which is needed for later mapping of RTT to
        % GetSecs time. In uncalibrated mode, this will just send out the
        % reset code.
        ResetRTT(handle);

        % Return handle:
        varargout{1} = handle;
%    catch
        % Close serial control link:
%        CloseDev(handle);
%    end
    return;
end

if strcmpi(cmd, 'GetBoxTimerSlope')
    % Close device:
    if nargin < 2
        error('You must provide the device "handle" for the box to compute slope for!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    % Perform calibrated basetimer query:
    [hosttime, devicetime] = queryBaseTimer(handle);
    
    lastBaseTimeQuery = ptb_cedrus_devices{handle}.lastBaseTimeQuery;
    
    baseToPtbSlope = (hosttime - lastBaseTimeQuery(1)) / (devicetime - lastBaseTimeQuery(2));
    
    % Store measured slope internally:
    ptb_cedrus_devices{handle}.baseToPtbSlope = baseToPtbSlope;

    % Return measured slope:
    varargout{1} = baseToPtbSlope;
    
    return;
end


if strcmpi(cmd, 'Close')
    % Close device:
    if nargin < 2
        error('You must provide the device "handle" for the box to close!');
    end

    % Retrieve handle and check if valid:
    handle = checkHandle(varargin{1});

    CloseDev(handle);
    
    return;
end

if strcmpi(cmd, 'CloseAll')
    % Close all open devices:

    if exist('ptb_cedrus_devices', 'var') 
        for handle = 1:length(ptb_cedrus_devices)
            if ~isempty(ptb_cedrus_devices{handle})
                CloseDev(handle);
            end
        end
        
        % All handles closed: Release the device array itself:
        clear ptb_cedrus_devices;
    end
    
    return;
end

error('Invalid subcommand given. Read the help.');
% ---- End of main routine ----


% ---- Start of internal helper functions ----

function TestThis(handle)
% Generic test blurb...
global ptb_cedrus_devices; %#ok<NUSED>

persistent testbyte;
if isempty(testbyte)
    testbyte = 0;
end

% Flush input buffer:
WaitSecs(0.2);
FlushEvents(handle);
WaitSecs(0.2);

% Set general mode for lines: "General purpose"
WriteDev(handle, 'a10');
WaitSecs(0.4);

% Set direction for lines: All output.
WriteDev(handle, ['a4' 0]);
WaitSecs(0.4);


% Query current state of outputs and inputs:
WriteDev(handle, 'ar');

% Wait for response:
inputLines = dec2bin(ReadDev(handle, 1)); %#ok<NASGU>

basetime = WaitSecs(0.5); %#ok<NASGU>

% Set all output lines low:
testbyte = mod(testbyte + 1, 256);
WriteDev(handle, ['ah' testbyte]);
%WriteDev(handle, ['ah' 255]);
WaitSecs(0.4);
WriteDev(handle, '_ah');

% Any activity, e.g., events???
%while BytesAvailable(handle) == 0
%    fprintf('Nothing yet at %f secs...\n', GetSecs - basetime);
%end

WaitSecs(0.1);
response = ReadDev(handle, 4); %#ok<NASGU>
%response2 = ReadDev(handle, 6)

% if length(response) == 6
%     [evt,CedrusStatus] = ExtractKeyPressData(handle,response)
% end

return;

function tReset = ResetRTT(handle)
% Try to reset the reaction time timer to zero within a small time
% window, so we can associate "time zero" of the RT timer with the
% current GetSecs() time. This way, the RT timer will encode elapsed
% time since that GetSecs basetime. We can then later on compute the
% time of a keypress simply as sum of the GetSecs baseline time and the
% reported event time (== value of RT timer at time of button press).
% Retry up to 100 times if reset doesn't occur within 10ms.
global ptb_cedrus_devices;

% Calibrated reset?
if ptb_cedrus_devices{handle}.baseToPtbSlope ~= 0
    % Calibrated reset:

    % Flush input buffer:
    WaitSecs(1);
    ClearQueues(handle);
    WaitSecs(1);

    % Switch to realtime priority if not already there:
    oldPriority=Priority;
    if oldPriority < MaxPriority('GetSecs')
        Priority(MaxPriority('GetSecs'));
    end
    
    % Get porthandle:
    blocking = 1;
    ntrials = 5;
    s = ptb_cedrus_devices{handle}.link;
    t = zeros(2,ntrials);
    
    % Perform up to ntrials trials:
    for ic=1:ntrials
        
        % Wait some random fraction of a millisecond. This will desync us
        % from the USB duty cycle and increase the chance of getting a very
        % small time window between scheduling, execution and acknowledge
        % of the send operation:
        WaitSecs(rand / 1000 + 1);
        
        % Take pre-Write timestamp: Sync command not emitted before that time:
        % Write sync command, wait 'blocking' for write completion, store
        % completion time in t(2,ic). Send RTT reset command code 'e5':
        [nw t(2,ic), errmsg, t(1,ic)] = IOPort('Write', s, 'e5', blocking);

        % We know that sync command emission has happened at some time
        % after t(1,ic) and before t(2,ic). This by design of the USB
        % standard, host controllers and operating system USB stack. This
        % is the only thing we can take for granted wrt. timing, so the
        % "time window" between those two timestamps is our window of
        % uncertainty about the real host time when sync started. However,
        % on a well working system without massive system overload one can
        % be reasonably confident that the real emission of the sync
        % command happened no more than 1 msec before t(2,ic). That is a
        % soft constraint however - useful for computing the final estimate
        % for hosttime, but nothing to be taken 100% for granted.
        
        if nw~=2
            % Send op failed!
            fprintf('CedrusResponseBox: RTTReset: Warning! Sync token send operation to box failed!\n');
            t(1,ic) = 0;
            t(2,ic) = inf;
            continue;
        end

        confidencewindow = t(2,ic) - t(1,ic);
        if confidencewindow < 0.001
            break;
        end

        % Next trial...
    end

    % Restore priority
    if Priority ~= oldPriority
        Priority(oldPriority);
    end
    
    % For each measurement, the time window t(2,ic)-t(1,ic) defines kind of
    % a confidence interval for the "real" host system time when the sync
    % command was emitted. The measurement with the smallest time window is
    % the most accurate one. Find it and use it:
    minwin = t(2,ic) - t(1,ic);
    
    % On OS/X or Linux we could easily do with 2 msecs, as a 1.2 msecs
    % minwin is basically never exceeded. On MS-Windows however, 2.x
    % durations are not uncommon, so we need to slack this to 3.
    if minwin > 0.002
        fprintf('CedrusResponseBox: RTTReset: Warning! Confidence interval for clock sync is %f msecs - More than 2 msecs!\n', minwin);
    end
    
    % If the 'minwin' window is smaller than 1 msec, we subtract 0.5 the
    % length of it from the t(2,idx) timestamp as best estimate for
    % hosttime -- reasonable assuming a uniform distribution in the
    % 'minwin' interval. If 'minwin' is more than 1 msecs, we assume it
    % happened 0.5 msecs before t(2,idx) -- taking advantage of the
    % soft-constraint that the real write usually happens within 1 msec of
    % t(2,idx) on a normally loaded and well working system:
    hosttime = t(2,ic) - (min(minwin, 0.001)/2);
    
    WaitSecs(0.2);
    FlushEvents(handle);
    WaitSecs(0.2);
    
    ptb_cedrus_devices{handle}.baseToPtbOffset = hosttime;
    ptb_cedrus_devices{handle}.rttresetdelay = minwin;
    tReset = minwin;
    
else
    % Only uncalibrated fast reset requested:

    % Reset offset fields to invalid values:
    ptb_cedrus_devices{handle}.baseToPtbOffset = 0;
    ptb_cedrus_devices{handle}.rttresetdelay = -1;

    % Send reaction time timer reset code:
    dx1 = GetSecs;
    WriteDev(handle, 'e5');
    dx2 = GetSecs;

    % Return estimated time of when reset probably roughly happened:
    tReset = (dx1+dx2)/2;
end

return;

function ClearQueues(handle)
    while BytesAvailable(handle)>0
        % Read and discard all bytes:
        ReadDev(handle, BytesAvailable(handle));
        WaitSecs(0.5);
    end
return;

function CedrusStatus=FlushEvents(handle)
% JCAR: I modified this because I don't want to just throw away information
% about key-presses. If I do, I lose track of what the current button
% status is. So, I will read the information and check whether any of it
% is key-presses. If so, I will use this information to update the current
% status.
global ptb_cedrus_devices;

CedrusStatus = ptb_cedrus_devices{handle}.CedrusStatus;

while BytesAvailable(handle)>0
    % Read 1 byte
    response=ReadDev(handle, 1);

    % See if this is "k", indicating that a key press is following
    if char(response)=='k'
        % Seems to be an event packet: Read remaining 5 more bytes
        last5=ReadDev(handle, 5);
        response(1:6) = [response(1) last5];
        [evt,CedrusStatus] = ExtractKeyPressData(handle,response);
    else
        fprintf('CedrusResponseBox:FlushEvents: Warning invalid value %s [%i] instead of "k" received!\n', char(response), response);
    end
end

return;

function label = findbuttonlabel(numbr, handle)
% The response box labels buttons by rather arbitrary numbers.
% I thought it might be helpful to have something more descriptive.
% THese descriptions assume the box is postioned with its cables/ports
% on the back edge furthest from the user.
global ptb_cedrus_devices;

switch (ptb_cedrus_devices{handle}.modelNo)
    case 530
        switch (numbr - 1)
            case 1
                label = 'top';
            case 6
                label = 'bottom';
            case 3
                label = 'left';
            case 5
                label = 'right';
            case 4
                label = 'middle';
            otherwise
                label = 'unknown';
        end
    case 730
        switch (numbr)
            case 2
                label = '1.Left';
            case 3
                label = '2.Left';
            case 4
                label = '3.Left';
            case 5
                label = '4.Left';
            case 6
                label = '5.Left';
            case 7
                label = '6.Left';
            case 8
                label = '7.Left';
            otherwise
                label = 'unknown';
        end
    % arrangement of key codes for Cedrus 830
    %   left    right   hand
    %   4  5    6  7
    % 8              1
    %       2  3
    case 830
        switch (numbr)
            case 1
                label = 'right.outer';
            case 2
                label = 'left.lower';
            case 3
                label = 'right.lower';
            case 4
                label = 'left.center';
            case 5
                label = 'left.inner';
            case 6
                label = 'right.inner';
            case 7
                label = 'right.center';
            case 8
                label = 'left.outer';
            otherwise
                label = 'unknown';
        end
    otherwise
        label = 'unknown';
end

return;

function [evt,CedrusStatus] = ExtractKeyPressData(handle,response)
%The XID device sends six bytes of information in the following format:
%<�k�><key info><RT>:
%
% The first parameter is simply the letter "k", lower case.
%
% The second parameter consists of one byte, divided into the following
% bits
%
% Bits 0-3 store the port number. For Lumina LP-400, the push buttons and
% scanner trigger are on port 0; the RJ45 I/O lines are on port 1.
% For SV?1, voice key is on port 2 and the RJ45 is on port 1 ? there is
% no port 0. For the RB-x30 response pads, the push buttons are on port 0
% and the RJ45 port is on port 1.
%
% Bit 4 stores an action flag. If set, the button has been pressed. If
% cleared, the button has been released.
%
% Bits 5-7 indicate which push button was pressed.
%
% The reaction time consists of four bytes and is the time elapsed since
% the Reaction Time timer was last reset. See description of command "e5".
%
% Information taken from http://www.cedrus.com/xid/protocols.htm

global ptb_cedrus_devices;

if length(response)~=6
    % Did not receive 6 bytes - This should not happen!
    error('In GetButtons: Received too short (or no) response packet from box!');
end

% According to cedrus, http://www.cedrus.com/xid/protocols.htm,
% The XID device sends six bytes of information in the following
% format: <�k�><key info><RT>:
% So the first byte is 107, ie the letter �k�, lower case

% Check byte 1 for correct value 'k':
if char(response(1))~='k'
    % Failed!
    error('Received invalid event packet [Not starting with a k] from box!');
end

% Extracts byte 2 to determine which button was pushed:
evt.raw = (response(2));

% Extract different bits into meaningful fields:
% According to Cedrus, the second parameter consists of one byte, divided into the following bits:
% Bits 0-3 store the port number.
% For Lumina LP-400, the push buttons and scanner trigger are on port 0;
% the RJ45 I/O lines are on port 1. For SV?1, voice key is on port 2 and
% the RJ45 is on port 1 � there is no port 0.
% For the RB-x30 response pads, the push buttons are on port 0 and the
% RJ45 port is on port 1.

% Port id: Bits 0-3
evt.port = bitand(evt.raw, 15);
% 15 is 1111, so this extracts the rightmost 4 bits from evt.raw, ie bits 0-3

% Button state: 1 = pressed, 0 = released. Bit 4
evt.action = bitand(bitshift(evt.raw, -4), 1);

% Button id: Which button? Bits 5-7
evt.button = bitshift(evt.raw, -5) + 1;
% This chops off the rightmost 5 bits, i.e. bits 0-4, leaving only bits
% 5-7

% Map to a more descriptive label: 
evt.buttonID = findbuttonlabel(evt.button, handle);

% Extracts bytes 3-6 and is the time elapsed in milliseconds since the
% Reaction Time timer was last reset.
%
% For more information about the use of XID timers refer to
% http://www.cedrus.com/xid/timing.htm
response = double(response);
evt.rawtime = 0.001 * (response(3)+(response(4)*256)+(response(5)*65536) +(response(6)*(65536*256)));

% Map rawtime to ptbtime if possible:
ptbTime = mapRTTimerToPTBTime(evt.rawtime, handle);

% Valid mapping? Assign if so. If mapping is impossible due to skipped
% timecalibration, we don't return the 'ptbtime' field. This way, usercode
% that relies on it without performing the mandatory calibration will die
% with a nice error message.
if ~isnan(ptbTime)
    evt.ptbtime = ptbTime;
end

% Try and keep track of which buttons are currently down and up, based on
% what bytes have been read in.
ptb_cedrus_devices{handle}.CedrusStatus(evt.port + 1, evt.button) = evt.action;

% CedrusStatus will tell you what buttons are currently up or down,
% based on the last time the device was read.
CedrusStatus = ptb_cedrus_devices{handle}.CedrusStatus;

return;

% Check if 'handle' is a valid handle into our struct array of devices,
% return it if it is valid, abort with error otherwise.
function retHandle = checkHandle(handle)
    global ptb_cedrus_devices;
    
    if handle > length(ptb_cedrus_devices) || isempty(ptb_cedrus_devices{handle})
        error('Invalid response box handle %i passed: No such response box device open!', handle);
    end
    retHandle = handle;
return;

% Helper function: Open serial connection:
function dev = OpenDev(port, baudrate)
global ptb_cedrus_drivertype;

% Test our default of type 2 -- Our own IOPort() driver:
if ptb_cedrus_drivertype == 2
    % Use IOPort:

    % Temporarily shut up the driver, so errors can be reasonably
    % handled:
    oldverb = IOPort('Verbosity', 0);

    % Open link:
     [dev.link, errmsg] = IOPort('OpenSerialPort', port, sprintf('BaudRate=%i Parity=None DataBits=8 StopBits=1 FlowControl=Hardware ReceiveTimeout=1 ', baudrate));

    IOPort('Verbosity', oldverb);

    % Success?
    if dev.link < 0
        % Nope. Do we know the cause?
        error(sprintf('Failed to open port %s for Cedrus response box via IOPort()! Reason: %s', port, errmsg)); %#ok<SPERR>
    end

    % Link is online.
    try
        
        % Clear all send and receive buffers and queues:
        IOPort('Purge', dev.link);

        % Assign output port, driverid and empty recvQueue:
        dev.driver = 2;
        dev.recvQueue = [];
        dev.port = port;

    catch
        error('Failed to open port %s for Cedrus response box via IOPort() driver.', port);
    end
    
    % Ready.
    return;
end

% Some non-standard driver: We support serial() on Windows and Linux,
% SerialComm on OS/X:
% Which OS?
if IsOSX
    % SerialComm:
    try
        % Open 'port' with 'baudrate' baud, no parity, 8 data bits, 1
        % stopbit.
        SerialComm('open', port, sprintf('%i,n,8,1', baudrate));

        % Disable handshaking 'n' == none:
        SerialComm('hshake', port, 'n');

        % Wait a bit...
        WaitSecs(0.5);

        % And flush all send- and receivebuffers:
        purgedata = SerialComm('read', port);

        if ~isempty(purgedata)
            fprintf('CedrusResponseBox: Open: Purged some trash data...\n');
        end

        % Assign and init stuff:
        dev.port = port;
        dev.link = port;
        dev.driver = 1;
        dev.recvQueue = [];
    catch
        error('Failed to open port %i on OS/X for Cedrus response box via SerialComm() driver.', port);
    end
else
    % Windows or Linux: Matlab supports serial() object in JVM mode:
    if ~psychusejava('desktop')
        error('You must run Matlab in JVM mode (JAVA enabled) for Cedrus response box to work!');
    end

    try
        % Ok, Matlab with JVM on Windows or Linux: Let's do it!
        dev.link = serial(port, 'BaudRate', baudrate, 'DataBits', 8, 'StopBits', 1,...
            'FlowControl', 'none', 'Parity', 'none', 'Terminator', 'CR', 'Timeout', 400,...
            'InputBufferSize', 16000);

        fopen(dev.link);
        dev.driver = 0;
        dev.port = port;
        dev.recvQueue = [];
    catch
        error('Failed to open port %s on Windows or Linux for Cedrus response box via Matlab serial() driver.', port);
    end
end

% Ready.
return;

function CloseDev(handle)
    global ptb_cedrus_devices;

    % Give device time to settle:
    WaitSecs(0.5);

    % Initiate a device reset:
    % WriteDev(handle, 'f7');

    % Give device time to settle after reset:
    WaitSecs(0.5);

    if ptb_cedrus_devices{handle}.driver == 0
        % Matlabs serial() driver:

        % Close serial control link:
        dev = ptb_cedrus_devices{handle};

        % Close data link:
        fclose(dev.link);

        % Delete serial control link object:
        delete(dev.link);
        clear dev.link;
    else
        if ptb_cedrus_devices{handle}.driver == 1
            % OS/X + Matlab + SerialComm driver:
            SerialComm('purge', ptb_cedrus_devices{handle}.link);
            SerialComm('close', ptb_cedrus_devices{handle}.link);
        end

        if ptb_cedrus_devices{handle}.driver == 2
            % IOPort driver:
            IOPort('Purge', ptb_cedrus_devices{handle}.link);
            IOPort('Close', ptb_cedrus_devices{handle}.link);
        end
    end
    
    % Clear out device struct:
    ptb_cedrus_devices{handle} = [];

return;

function nrAvail = BytesAvailable(handle)
    global ptb_cedrus_devices;

    if ptb_cedrus_devices{handle}.driver == 0
        % Matlabs serial() driver:
        
        % Readout BytesAvailable subfield of device link object:
        nrAvail = ptb_cedrus_devices{handle}.link.BytesAvailable;
    else
        if ptb_cedrus_devices{handle}.driver == 1
            % OS/X + Matlab + SerialComm driver:
            
            % All reads are non-blocking and there isn't any BytesAvailable
            % command. We fetch all data that's currently available via
            % non-blocking read and attach it to our own queue, then return
            % the total number of bytes in the queue:
            data = transpose(SerialComm('read', ptb_cedrus_devices{handle}.link));
            ptb_cedrus_devices{handle}.recvQueue = [ptb_cedrus_devices{handle}.recvQueue data];
            nrAvail = length(ptb_cedrus_devices{handle}.recvQueue);
        end

        if ptb_cedrus_devices{handle}.driver == 2
            % IOPort driver:            
            nrAvail = IOPort('BytesAvailable', ptb_cedrus_devices{handle}.link);
        end
    end
    
    if nrAvail > 0
        % Store timestamp when queue was not empty:
        ptb_cedrus_devices{handle}.lastTimeQueueNonEmpty = GetSecs;
    end
    
return;

function data = ReadDev(handle, nwanted)
    global ptb_cedrus_devices;

    if ptb_cedrus_devices{handle}.driver == 0
        % Matlabs serial() driver:
        
        % Read via fread 'nwanted' bytes from link. Block until we get the
        % wanted 'nwanted' bytes or until timeout / error:
        data = transpose(fread(ptb_cedrus_devices{handle}.link, nwanted));
    else
        if ptb_cedrus_devices{handle}.driver == 1
            % OS/X + Matlab + SerialComm driver:
            
            % Call BytesAvailable to trigger read-in of data from serial
            % port to our internal queue and to update the available stats,
            % until at least the 'nwanted' bytes are available, or until
            % the read operation times out after 2 seconds:
            currtime = GetSecs;
            timeout = currtime + 2;
            while (BytesAvailable(handle) < nwanted) && (currtime < timeout)
                % We are on OS/X, so waiting for 1 msec should suffice, no
                % need to wait 4 msecs as on that other deficient OS:
                currtime = WaitSecs(0.001);
            end;
            
            if currtime >= timeout
                fprintf('Timed out: nwanted = %i, got %i bytes: %s\n', nwanted, BytesAvailable(handle), char(ptb_cedrus_devices{handle}.recvQueue));
                fprintf('Read operation on response box timed out after 2 secs!\n');
                data = [];
                return;
            end
            
            % Have at least the nwanted bytes, so fetch the first nwanted
            % bytes from queue:
            data = ptb_cedrus_devices{handle}.recvQueue(1:nwanted);
            
            % Dequeue them from queue:
            if length(ptb_cedrus_devices{handle}.recvQueue) > nwanted
                % Keep tail of queue:
                ptb_cedrus_devices{handle}.recvQueue = ptb_cedrus_devices{handle}.recvQueue(nwanted+1:end);
            else
                % Nothing more in queue: Delete it.
                ptb_cedrus_devices{handle}.recvQueue = [];
            end            
        end

        if ptb_cedrus_devices{handle}.driver == 2
            % IOPort driver: Returns all data as data type double:
            % fprintf('In read....\n');
            [data, when, errmsg] = IOPort('Read', ptb_cedrus_devices{handle}.link, 1, nwanted);
            if length(data) < nwanted
                fprintf('Timed out: nwanted = %i, got %i bytes: %s\n', nwanted, length(data), char(data));
                fprintf('Read operation on response box timed out after 1 secs! errmsg = %s\n', errmsg);
                data = [];
                return;
            end
        end
    end
return;

function WriteDev(handle, data)
    global ptb_cedrus_devices;

    if ptb_cedrus_devices{handle}.driver == 0
        % Matlabs serial() driver:
        
        % Write data via fwrite: We provide our own '%s' formatting string
        % to make sure that data is passed as-is, without any terminators
        % (CR or LF or CR+LF) attached. This will block until send
        % completion:
        fwrite(ptb_cedrus_devices{handle}.link, char(data));
    else
        if ptb_cedrus_devices{handle}.driver == 1
            % OS/X + Matlab + SerialComm driver:
            
            % Write data - without terminator - via SerialComm:
            SerialComm('write', ptb_cedrus_devices{handle}.link, double(data));
        end

        if ptb_cedrus_devices{handle}.driver == 2
            % IOPort driver:
            
            % Write data - without terminator:
            % fprintf('In write....\n');
            IOPort('Write', ptb_cedrus_devices{handle}.link, char(data), 1);
        end
    end
return;

function roundtrip = RoundTripTestDev(handle)
    global ptb_cedrus_devices;

    % Flush the queue:
    FlushEvents(handle);

    % Perform 100 measurement trials:
    roundtrip = zeros(1,100);
    for i=0:100
        % Wait a bit between each trial:
        WaitSecs(0.100);

        % Send 'e4' code to initiate procedure:
        WriteDev(handle, 'e4');

        % Wait for receive completion:
        while BytesAvailable(handle) < 1
        end;

        % Send echo, optimistically assuming we received a 'X':
        WriteDev(handle, 'X');

        % Get the really received byte and check:
        if char(ReadDev(handle, 1))~='X'
            error('Roundtrip test did not receive "X" char as expected!');
        end

        % Wait for receipt of timestamp:
        while BytesAvailable(handle) < 4
        end;

        response = ReadDev(handle, 4);

        if response(1)~='P' || response(2)~='T'
            error('Roundtrip test did not receive "PT" marker as expected!');
        end

        response = double(response);
        
        % We throw away the first trial:
        if i > 0
            roundtrip(i) = 0.001 * (response(3) + 256 * response(4));
        end
    end

    % Store median and stddev of roundtrip time in device struct:
    ptb_cedrus_devices{handle}.roundtriptime = median(roundtrip);
    ptb_cedrus_devices{handle}.roundtripstddev = std(roundtrip);

return;

function [hosttime, devicetime, minwin] = queryBaseTimer(handle)
    global ptb_cedrus_devices;

    % Flush input buffer:
    WaitSecs(1);
    ClearQueues(handle);
    WaitSecs(1);

    % Switch to realtime priority if not already there:
    oldPriority=Priority;
    if oldPriority < MaxPriority('GetSecs')
        Priority(MaxPriority('GetSecs'));
    end
    
    % Get porthandle:
    blocking = 1;
    ntrials = 5;
    s = ptb_cedrus_devices{handle}.link;
    t = zeros(2,ntrials);
    
    % Perform up to ntrials trials:
    for ic=1:ntrials
        
        % Wait some random fraction of a millisecond. This will desync us
        % from the USB duty cycle and increase the chance of getting a very
        % small time window between scheduling, execution and acknowledge
        % of the send operation:
        WaitSecs(rand / 1000 + 1);
        
        % Take pre-Write timestamp: Sync command not emitted before that time:
        % Write sync command, wait 'blocking' for write completion, store
        % completion time in t(2,ic). Send basetimer query command code 'e3':
        [nw t(2,ic), errmsg, t(1,ic)] = IOPort('Write', s, 'e3', blocking);

        % Wait for response from box, receive packet, 
        % then parse into raw timer value (in seconds):
        devicetime = receiveAndParseTimePacket(handle);
        
        % We know that query command emission has happened at some time
        % after t(1,ic) and before t(2,ic). This by design of the USB
        % standard, host controllers and operating system USB stack. This
        % is the only thing we can take for granted wrt. timing, so the
        % "time window" between those two timestamps is our window of
        % uncertainty about the real host time when sync started. However,
        % on a well working system without massive system overload one can
        % be reasonably confident that the real emission of the sync
        % command happened no more than 1 msec before t(2,ic). That is a
        % soft constraint however - useful for computing the final estimate
        % for hosttime, but nothing to be taken 100% for granted.
        if nw~=2
            % Send op failed!
            fprintf('CedrusResponseBox: queryBaseTimer: Warning! Query token send operation to box failed!\n');
            t(1,ic) = 0;
            t(2,ic) = inf;
            continue;
        end

        confidencewindow = t(2,ic) - t(1,ic);
        if confidencewindow < 0.001
            break;
        end

        % Next trial...
    end

    % Restore priority
    if Priority ~= oldPriority
        Priority(oldPriority);
    end
    
    % For each measurement, the time window t(2,ic)-t(1,ic) defines kind of
    % a confidence interval for the "real" host system time when the sync
    % command was emitted. The measurement with the smallest time window is
    % the most accurate one. Find it and use it:
    minwin = t(2,ic) - t(1,ic);
    
    % On OS/X or Linux we could easily do with 2 msecs, as a 1.2 msecs
    % minwin is basically never exceeded. On MS-Windows however, 2.x
    % durations are not uncommon, so we need to slack this to 3.
    if minwin > 0.002
        fprintf('CedrusResponseBox: queryBaseTimer: Warning! Confidence interval for clock sync is %f msecs - More than 2 msecs!\n', minwin);
    end
    
    % If the 'minwin' window is smaller than 1 msec, we subtract 0.5 the
    % length of it from the t(2,idx) timestamp as best estimate for
    % hosttime -- reasonable assuming a uniform distribution in the
    % 'minwin' interval. If 'minwin' is more than 1 msecs, we assume it
    % happened 0.5 msecs before t(2,idx) -- taking advantage of the
    % soft-constraint that the real write usually happens within 1 msec of
    % t(2,idx) on a normally loaded and well working system:
    hosttime = t(2,ic) - (min(minwin, 0.001)/2);
    
    WaitSecs(0.2);
    FlushEvents(handle);
    WaitSecs(0.2);
    
return;

% Reads raw basetimer response packet from box, converted to seconds, but
% not corrected for receive latency etc. Query command must have been sent
% by calling code!
function rawBaseTime = receiveAndParseTimePacket(handle)

    % Read all 6 bytes of basetimer response packet from box:
    response = ReadDev(handle, 6);

    if length(response)~=6
        % Did not receive 6 bytes - This should not happen!
        error('In receiveAndParseTimePacket: Received too short (or no) response packet from box!');
    end

    % Check bytes 1:2 for correct values 'e3':
    if char(response(1))~='e' || char(response(2))~='3'
        % Failed!
        error('In receiveAndParseTimePacket: Received invalid response packet [Not starting with "e3"] from box!');
    end

    % Extracts bytes 3-6 and is the time elapsed in milliseconds since the
    % base timer was last reset.
    %
    % For more information about the use of XID timers refer to
    % http://www.cedrus.com/xid/timing.htm
    % Conver to seconds:
    response = double(response);
    rawBaseTime = 0.001 * (response(3)+(response(4)*256)+(response(5)*65536) +(response(6)*(65536*256)));

return;

function ptbTime = mapRTTimerToPTBTime(rtt, handle)
    global ptb_cedrus_devices;

    if ptb_cedrus_devices{handle}.baseToPtbOffset ~= 0
        % rtt is the parsed timevalue (already mapped from msecs to seconds),
        % as received in a event packet from the box. We map it to ptbTime by
        % adding the offset between GetSecs time and device RTT time, as
        % estimated by last calibrated RTTReset():
        ptbTime = ptb_cedrus_devices{handle}.baseToPtbOffset + rtt;
    else
        % Missing clock sync. Return "invalid" result:
        ptbTime = nan;
    end
    
return;
