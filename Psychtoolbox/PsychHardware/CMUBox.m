function varargout = CMUBox(cmd, handle, varargin)
% CMUBox - Access CMU response button box or PST serial response button box as well as fORP and Bitwhacker devices.
% 
%
% This allows to query button response boxes of type CMU (Carnegie Mellon
% University box) and PST (E-Prime response box). It also allows to use a
% UBW32/Bitwhacker device to be used as a response box if the device is
% loaded with the StickOS firmware from http://cpustick.com. It also allows
% to use the Curdes fORP devices if connected via serial port. And it offers
% a simple way to access the RTBox in E-Prime mode.
%
% Commands and their syntax:
% --------------------------
%
% handle = CMUBox('Open', boxtype [, portName] [, options] [, debounceSecs=0.030] [, isInverted]);
% - Open response box connected to serial port 'portName', or the first
% serial port found, if 'portName' is omitted. Initialize it, return a
% 'handle' to it. You'll have to pass 'handle' to all following functions
% to access the box.
%
% If your system has multiple devices connected to multiple serial ports
% then you should explicitely specify the 'portName', otherwise the driver
% may connect to the wrong port and choke!
%
% The optional string parameter 'options' allows to tweak the behaviour of
% the driver for certain configurations. It supports the following options:
%
% 'ftdi' - Tells the driver that it is connecting to a Serial-over-USB port
% and that the converter/driver is from FTDI Inc., or a compatible device.
% Allows for certain optimizations in timing accuracy.
%
% 'norelease' - Tells the driver it shouldn't report button release / TTL
% transition to low events, but only button presses and TTL onsets.
%
%
% The optional parameter 'debounceSecs' sets the debounce interval for
% button debouncing if the Bitwhacker is used as response box. After a
% button press, the device will wait for all the buttons to be released for
% at least 'debounceSecs' seconds. Only then will it accept new button
% presses. The default setting is 30 msecs if this option is omitted.
%
%
% The optional parameter vector 'isInverted' defines whether a TTL input
% signal level of logic low (0) corresponds to a button press or a button
% release. By default, a level of logic low is detected as button press.The
% vector has nine elements, one for each corresponding input line. If an
% element is 1, then a logic level of low is registered as a button press.
% A value of 0 means that a logic level of high is registered as a button
% press. E.g., isInverted = [ 0 0 0 0 1 1 1 1 1 ] would register a logic
% level of low as a button press on the first four inputs, whereas a level
% of high is required to detect a button press on the last five inputs. 
% 
%
% The mandatory parameter 'boxtype' is a name string defining the
% type/model of box to connnect to. Supported settings are:
%
% 'bitwhacker' - Connect to UBW32/Bitwhacker with StickOS Firmware,
% reconfigure it to impersonate a response box.
%
% 'cmu' - Connect to CMU serial port response button box, assuming 19.2
% KBaud datarate, 8 databits, 1 stop bit, odd parity.
%
% 'pst' - Connect to PST serial port response button box in E-Prime
% configuration with 800 streaming samples per second at 19200 Baud, 8
% databits, 1 stopbit, no parity. Please note this is our expectation: If
% the box is configured for something else than 19200 Baud, we will fail or
% hang. If 1600 samples/sec are configured instead of 800 sampes/sec, we
% won't care, but timing might be less accurate due to the too high USB load,
% unless we're running on a real serial port which shouldn't have a
% problem with 1600 samples/sec. Serial-over-USB however will likely choke!
% Verify that the DIP switches or jumpers inside the box are properly set
% up.
%
% 'forpserial-1' - Connect to a fORP interface unit type FIU-005 on the
% serial port, configured with a program switch setting of 1. This is
% exactly the same as a PST box running with 800 Samples/Sec streaming
% rate.
%
% 'forpserial-0', 'forpserial-2', 'forpserial-4', 'forpserial-6'
% Connect to a fORP interface unit type FIU-005 on the serial port,
% configured with a program switch setting of 0, 2, 4 or 6. In this mode,
% similar to the Bitwhacker, the device will only send a single byte if the
% status of the buttons or triggers changes. This is more efficient than
% mode 'forpserial-1', so one of these settings is recommended.
%
% Modes 0, 4 and 6 will only report button presses or trigger reception,
% but not button releases, whereas mode 2 will report any status change,
% ie., it will also report and timestamp button releases. The mapping of
% the returned status value to corresponding button / trigger states is
% dependent on the selected mode. Modes 0, 4 and 6 report ASCII codes that
% identify the pressed button, whereas mode 2 returns a status byte where
% each bit encodes the current (updated) status of a single button, similar
% to the PST and CMU response boxes. For the specifics, see the fORP manual
% at http://www.curdes.com/ForpUserGuide.html
%
% 'rtbox' - Connect to a RTBox in "simple mode", the mode it uses after
% powerup. The Box will send a byte of data for each event, each bit encoding
% the new status of one button or other input.
%
% 'lumina' - Connect to a Cedrus Lumina response box.
%
%
% CMUBox('Close', handle);
% - Close connection to response box 'handle'. The 'handle' is invalid
% thereafter.
%
%
% evt = CMUBox('GetEvent', handle [, waitForEvent=0]);
% - Retrieve next queued event received from the box in the struct 'evt'.
% If no new events are available and the optional 'waitForEvent' is set to
% 1, then the function will wait until at least one valid event becomes
% available and return that event. Otherwise it will return an empty struct,
% ie., evt = [] to signal that no new events are available.
%
% The following subfields are available in 'evt' if 'evt' is non-empty:
%
% evt.state = 1 Byte value which encodes the new status of the response
% buttons and input lines of the box. A 1-bit means button pressed/signal
% active. A 0-bit means button released/signal inactive. See the
% documentation of your box for meaning of the single bits. If you use the
% Bitwhacker then evt.state does directly encode the button number of a
% pressed button: 0 == All buttons released and all TTL inputs low.
% Values 1-7 correspond to a low->high transition of TTL input pins A1 - A7.
% Values 8 and 9 correspond to a button press of onboard buttons "USER" or
% "PRG". Bitwhacker can only report one active button or TTL line at a time.
%
% The fORP device will encode evt.state differently depending on selected
% mode, either like a Bitwhacker device or like a CMU/PST response box.
%
% evt.time  = Psychtoolbox GetSecs() timestamp of the time when the event
% was received from the box. The accuracy depends on the properties of your
% serial port device and system load. For a native serial port and a
% normally loaded system, you can expect about 1-2 msec delay. For a
% Serial-over-USB port, it depends on the converter type and driver
% settings but can be 2-3 msecs at best.
%
% evt.trouble = If zero, then evt is probably valid and good. If non-zero,
% then the timestamp is likely screwed and useless, as are probably all
% following timestamps, at least for CMU/PST boxes and fORP's in mode 1!
% For the Bitwhacker or fORP's in modes 0,2,4 or 6, only the current 'evt'
% will be invalid, but later events will recover and thereby be unaffected.
%
%
% status = CMUBox('Status', handle);
% - Retrieve internal status of response box 'handle' as a struct.
%
%

% How to setup latency timer and eventChar under Linux for FTDI
% USB-Serial converters?
%
% Example for Ubuntu Linux 7.1 with first serial converter device ttyUSB0:
%
% su + Superuser password.
% cd /sys/bus/usb-serial/devices/ttyUSB0
%
% Set event character:
%
% echo X > event_char
% with X being 0 to disable event character, or 256 + ASCII code of
% event character, e.g., echo 256 > event_char for ASCII code 0, or
% echo 266 > event_char for ASCII code 10, because 266 = 256 + 10.
%
% Set latency timer to 1 millisecond: (Default is 16 msecs)
% echo 1 > latency_timer
%
% Query current latency timer setting:
% cat latency_timer
%
% History:
% 9.08.2009  mk  Written. Initial prototype.
% 3.03.2010  mk  Update to properly use Bitwhacker as response box as well
%                as the Curdes fORP devices when connected via serial port.
% 13.02.2012 mk  Update with support for Cedrus Lumina, based on
%                information provided by some user named "Nick".

% Cell array of structs for our boxes: One cell for each open box.
persistent boxes;

if nargin < 1
    help CMUBox;
    return;
end

if nargin < 2
    handle = [];
end

if strcmpi(cmd, 'GetEvent')
    if isempty(handle)
        error('CMUBox: GetEvent: No "handle" for box to close provided!');
    end

    if length(handle) ~= 1
        error('CMUBox: GetEvent: Passed argument is not a valid response box "handle"!');
    end
    
    if length(boxes) < handle
        error('CMUBox: GetEvent: Passed argument is not a valid response box "handle"!');
    end

    % Retrieve box struct:
    box = boxes{handle};
    
    if isempty(box.useBitwhacker)
        error('CMUBox: GetEvent: Passed argument is not a valid response box "handle"!');
    end
    
    
    if length(varargin) < 1
        waitEvent = [];
    else
        waitEvent = varargin{1};
    end
    
    if isempty(waitEvent)
        waitEvent = 0;
    end
    
    % Start with empty event 'evt':
    evt = [];
    
    % Reset timestamp trouble flag:
    tTrouble = 0;
    skipped = 0;
    
    % Repeat until forceful abortion via break as long as either waiting
    % for new events is requested, or - in non-blocking mode - new status
    % bytes from the box are available to parse:
    while (waitEvent > 0) || (IOPort('BytesAvailable', box.port) >= 9)

        % Wait blocking for at least one status packet of 9 bytes from box:
        [inpkt, t, err] = IOPort('Read', box.port, 1, 9);

        % Error condition?
        if ~isempty(err)
            % Game over!
            error('CMUBox: GetEvent: I/O ERROR!! System says: %s\n', err);
        end

        % No data returned due to read timeout? Just repeat loop.
        if isempty(inpkt)
            continue;
        end
        
        % No valid data returned due to read timeout? Just repeat loop.
        if ~any(inpkt)
            continue;
        end
        
        % Box status byte:
        data = inpkt(1);
        
        % Serial number of byte in stream:
        serNumber = inpkt(2) * 256^0 + inpkt(3) * 256^1 + inpkt(4) * 256^2 + inpkt(5) * 256^3;
        refTime   = (serNumber * box.dt) + box.baseTime;
        
        % Time delta between consecutive stream scans in seconds:
        % Effectively the sampling interval and therefore uncertainty of
        % the time measurement: If this deviates significantly from the
        % expected per-byte transmission interval then the timestamps are
        % not to be trusted!
        box.deltaScan = (inpkt(6) * 256^0 + inpkt(7) * 256^1 + inpkt(8) * 256^2 + inpkt(9) * 256^3) / 1e6;

        % Special case for Bitwhacker emulation: Filter out codes 10 and
        % 13, as well as 0 as they're an artifact of the emulation:
        if box.useBitwhacker && ismember(data, [0, 10, 13])
            continue;
        end
        
        % Timestamps at least 0.5 msecs apart and no more than 2 msecs
        % apart? This window should be sufficient for the CMU and PST box
        % in all streaming modes:
        if (t - box.oldTime < 0.0005) || (box.deltaScan < 0.0005) || ((box.deltaScan > 0.002) && (box.Streaming > 0))
            % Too close to each other! Timestamp is not reliable!
            tTrouble = 1;
            fprintf('CMUBox: GetEvent: Timestamp trouble!! Delta %f msecs, ScanInterval %f msecs.\n', 1000 * (t - box.oldTime), 1000 * box.deltaScan); 
        end

        % Keep track of last events timestamp:
        box.oldTime = t;
        
        % Update last status state of box:
        oldState = box.oldState;
        box.oldState = data;

        % Update trouble counter:
        box.tTrouble = box.tTrouble + tTrouble;

        % Box status changed since last query?
        if (data ~= oldState) 
            % Yes. We have a new event. Store it and break out of loop:
            
            % USB-Serial converter type?
            if box.ftdiusbserialtype > 0
                % Special timestamp correction for some transitions possible:
                if oldState == 0
                    % Zero to non-zero transition: As zeros are reported with 1
                    % msec internal latency, and the latency timer is reset
                    % at each reported zero, the first non-zero will be
                    % delayed by the latency of the FTDI latencytimer,
                    % i.e., with 2 msecs delay. Therefore we can subtract 2
                    % msecs to compensate for this case:
                    t = t - 0.002;
                else
                    if data == 0
                        % Non-zero to zero transition: Transmitted with
                        % only 1 msec delay due to FTDI special event
                        % character set to trigger on zeros, as requested
                        % by setting the 'Terminator' to 0. We subtract 1
                        % msecs to take this into account:
                        t = t - 0.001;
                    end
                end
            end

            % For Bitwhacker, remap ASCII character codes for symbols "0"
            % to "9" back to values 0 to 9. A 0 means "No button pressed",
            % whereas 1-9 mean button 1 to 9 pressed:
            if box.useBitwhacker
                data = data - 48;
            end
            
            % Should we discard "release" events and this is one?
            if (box.norelease == 0) || (data ~= 0)
                % Nope. Either no release event or don't discard'em.
                % Build 'evt' struct and return:                

                if ~box.Streaming
                    % streamTime is meaningless in non-streaming mode:
                    refTime = 0;
                end
                
                evt.time = t;
                evt.streamTime = refTime;
                evt.state = data;
                evt.trouble = tTrouble;
                evt.deltaScan = box.deltaScan;
                break;
            else
                % This is a release event that shall be discarded:
                evt = [];
            end
        end
        
        % Repeat scanloop until we find a new event or no new data
        % available from box in non-blocking mode...
        skipped = skipped + 1;
    end

    box.skipped = skipped;
    
    % Update box status:
    boxes{handle} = box;
    
    % Return evt if any:
    varargout{1} = evt;
    
    return;
end

% Open connection to box:
if strcmpi(cmd, 'Open')
    if isempty(handle)
        error('CMUBox: Open: Mandatory "boxtype" parameter missing!');
    else
        boxtype = handle;
    end

    %IOPort('Verbosity', 9);

    switch(lower(boxtype))
        case {'bitwhacker'},
            % No special options for UBW32/Bitwhacker:
            pString = '';
            box.useBitwhacker = 1;
            box.Streaming = 0;
            box.type = 1;
            fprintf('CMUBox: Using Bitwhacker/StickOS emulated box!\n');
            
        case {'forpserial-0', 'forpserial-4', 'forpserial-6'},
            % % BaudRate is 19.2 KiloBaud, 8-N-1 config without flow
            % control:
            pString = 'BaudRate=19200';
            box.useBitwhacker = 0;
            box.Streaming = 0;
            box.type = 1;
            fprintf('CMUBox: Using fORP interface program "%s" as serial response button box!\n', lower(boxtype));

        case {'forpserial-2'},
            % % BaudRate is 57.6 KiloBaud, 8-N-1 config without flow
            % control:
            pString = 'BaudRate=57600';
            box.useBitwhacker = 0;
            box.Streaming = 0;
            box.type = 1;
            fprintf('CMUBox: Using fORP interface program 2 as serial response button box!\n');
            
        case {'rtbox'},
            % % BaudRate is 115.2 KiloBaud, 8-N-1 config without flow
            % control:
            pString = 'BaudRate=115200';
            box.useBitwhacker = 0;
            box.Streaming = 0;
            box.type = 1;
            fprintf('CMUBox: Using RTBox as serial response button box!\n');
            
        case {'cmu'},
            % BaudRate is 19.2 KiloBaud, 8-Odd-1 config without flow control:
            % This is valid for at least hardware type mark 3, and mark 4,
            % Firmware V.1.x.
            pString = 'BaudRate=19200 Parity=Odd ReceiveTimeout=10.0';
            box.useBitwhacker = 0;
            box.Streaming = 1;
            box.type = 2;
            fprintf('CMUBox: Using CMU serial response button box!\n');
            
        case {'pst'},
            % BaudRate is 19.2 KiloBaud, 8-N-1 config without flow control:
            pString = 'BaudRate=19200 ReceiveTimeout=10.0';
            box.useBitwhacker = 0;
            box.Streaming = 1;
            box.type = 3;
            fprintf('CMUBox: Using PST serial response button box!\n');
            
        case {'forpserial-1'},
            % BaudRate is 19.2 KiloBaud, 8-N-1 config without flow control:
            pString = 'BaudRate=19200 ReceiveTimeout=10.0';
            box.useBitwhacker = 0;
            box.Streaming = 1;
            box.type = 3;
            fprintf('CMUBox: Using fORP interface program 1 as serial response button box!\n');

        case {'lumina'},
            % BaudRate is 9.6 KiloBaud, 8-E-1 config without flow control:
            pString = 'BaudRate=9600 Parity=Even';
            box.useBitwhacker = 0;
            box.Streaming = 0;
            box.type = 1;
            fprintf('CMUBox: Using Cedrus Lumina box!\n');

        otherwise,
            error('CMUBox: Open: Unknown "boxtype" specified! Typo?');
    end
    
    if length(varargin) < 1
        portName = [];
    else
        portName = varargin{1};
    end
    
    if isempty(portName)
        % No portName provided: Try to auto-detect:
        portName = FindSerialPort([], 1);
    end
    
    if length(varargin) >= 2
        specialOptions = varargin{2};
    else
        specialOptions = '';
    end
    
    if ~isempty(findstr(specialOptions, 'ftdi'))
        % USB-Serial port via FTDI converter chip:
        box.ftdiusbserialtype = 1;
    else
        % Default to non-FTDI:
        box.ftdiusbserialtype = 0;
    end

    if ~isempty(findstr(specialOptions, 'norelease'))
        % Don't report button release events:
        box.norelease = 1;
    else
        % Default is to report release events if supported:
        box.norelease = 0;
    end

    if IsLinux
        % Switch serial port to low-latency mode on Linux:
        % This is extra-paranoid, as our udev rules and the
        % 'ReceiveLatency=0.0001' setting below does already the same.
        system(sprintf('setserial %s low_latency', port));
        WaitSecs('YieldSecs', 0.100);
    end

    % Try to open connection: Allocate an input buffer of a size of
    % 1600 * 9 * 3600 = 51840000 Bytes. This is sufficient for 1 hour of
    % uninterrupted box operation without ever reading out events from the
    % queue at highest box operating speed.
    % On MS-Windows, we set the driver receivebuffer seize to 32 kB.
    %
    % The Terminator=0 setting is interesting for USB-Serial converters of
    % FTDI Inc. It will configure the converter to transmit zero bytes,
    % ie., "button resting position" immediately at each USB work-cycle
    % with at most 1 msec converter-internal latency. Other states, ie.,
    % "at least one button pressed", will be reported with the latency set
    % in the latency timer - at least 2 msecs delay. This however means
    % that a zero -> non-zero transition is known to have 2 msecs delay and
    % we can account for that, and that a non-zero -> zero transition has 1
    % msec delay and we can account for that. Net result is that we get a
    % low-latency for button transitions, as long as they are from none to
    % some and some to none.
    box.portName = portName;
    box.port = IOPort('OpenSerialPort', portName, ['InputBufferSize=51840000 HardwareBufferSizes=32768,32768 Terminator=0 ReceiveLatency=0.0001 ' pString]);
    
    
    % Is this a testrun with an emulated CMU/PST box by use of the
    % UBW32-Bitwhacker device with StickOS?
    if box.useBitwhacker
        % Upload StickOS Basic program to turn the BitWhacker into a CMU
        % box look-alike:

        % Upload simple program:
        % ----------------------

        % Disable input echoing from StickOS:
        Command(box.port, 'echo off');
        % Delete current program, upload a fresh one:
        Command(box.port, 'new');

        if 1
            if length(varargin) >= 3
                deadTimeSecs = varargin{3};
            else
                deadTimeSecs = [];
            end

            if isempty(deadTimeSecs)
                deadTimeSecs = 0.030;
            end

            deadTimeSecs = max(ceil(deadTimeSecs * 1000), 0);
            fprintf('Bitwhacker: Debounce: Will ignore new button presses for %i msecs after previous button release.\n', deadTimeSecs);
            
            if length(varargin) < 4
                isInverted = [];
            else
                isInverted = varargin{4};
            end
            
            if isempty(isInverted)
                isInverted = [1 1 1 1 1 1 1 1 1];
            else
                if length(isInverted) ~= 9 || ~isnumeric(isInverted)
                    error('Parameter "isInverted" must be a 9-element vector of zeros and ones!');
                end
            end
            
            for i=1:9
                if isInverted(i)
                    inv{i} = 'inverted '; %#ok<AGROW>
                else
                    inv{i} = ''; %#ok<AGROW>
                end
            end
            
            % Infinite while-loop, runs until program termination:
            Command(box.port, '70  dim led0 as pin re3 for digital output inverted');
            Command(box.port, '80  dim led1 as pin re2 for digital output inverted');
            Command(box.port, '90  dim busyled as pin re1 for digital output inverted');
            Command(box.port, ['100 dim ttl1 as pin ra1 for digital input ' inv{1}]);
            Command(box.port, ['101 dim ttl2 as pin ra2 for digital input ' inv{2}]);
            Command(box.port, ['102 dim ttl3 as pin ra3 for digital input ' inv{3}]);
            Command(box.port, ['103 dim ttl4 as pin ra4 for digital input ' inv{4}]);
            Command(box.port, ['104 dim ttl5 as pin ra5 for digital input ' inv{5}]);
            Command(box.port, ['105 dim ttl6 as pin ra6 for digital input ' inv{6}]);
            Command(box.port, ['106 dim ttl7 as pin ra7 for digital input ' inv{7}]);
            Command(box.port, ['107 dim but1 as pin re6 for digital input ' inv{8}]);
            Command(box.port, ['108 dim but2 as pin re7 for digital input ' inv{9}]);
            Command(box.port, '110 dim deadline');
            Command(box.port, '120 dim deadtime');
            Command(box.port, '130 dim ttlsum');
            Command(box.port, sprintf('140 let deadtime = %i', deadTimeSecs));
            Command(box.port, '160 print "0"');
            Command(box.port, '170 let led0 = 0');            
            Command(box.port, '180 let led1 = 0');            
            Command(box.port, '190 let busyled = 0');            
            Command(box.port, '200 while 1 do');
            Command(box.port, '210   let ttlsum = ttl1 * 1 + ttl2 * 2 + ttl3 * 3 + ttl4 * 4');
            Command(box.port, '220   let ttlsum = ttlsum + ttl5 * 5 + ttl6 * 6 + ttl7 * 7');
            Command(box.port, '230   let ttlsum = ttlsum + but1 * 8 + but2 * 9');
            Command(box.port, '300   if ttlsum > 0 then');
            Command(box.port, '400     print ttlsum');
            Command(box.port, '405     let busyled = 1');
            Command(box.port, '407     let led0 = but2 + ttl1 + ttl3 + ttl5 + ttl7');
            Command(box.port, '408     let led1 = but1 + ttl2 + ttl4 + ttl6');
            Command(box.port, '410     let deadline = msecs + deadtime');
            Command(box.port, '420     while msecs < deadline do');
            Command(box.port, '422       let ttlsum = ttl1 * 1 + ttl2 * 2 + ttl3 * 3 + ttl4 * 4');
            Command(box.port, '423       let ttlsum = ttlsum + ttl5 * 5 + ttl6 * 6 + ttl7 * 7');
            Command(box.port, '424       let ttlsum = ttlsum + but1 * 8 + but2 * 9');
            Command(box.port, '427       let led0 = but2 + ttl1 + ttl3 + ttl5 + ttl7');
            Command(box.port, '428       let led1 = but1 + ttl2 + ttl4 + ttl6');
            Command(box.port, '430       if ttlsum > 0 then');
            Command(box.port, '440         let deadline = msecs + deadtime');
            Command(box.port, '450       endif');
            Command(box.port, '460     endwhile');
            Command(box.port, '470     print "0"');
            Command(box.port, '480     let busyled = 0');
            Command(box.port, '500   endif');
            Command(box.port, '600 endwhile');
        else
            % Test code for MK's internal debugging:
            
            % Infinite while-loop, runs until program termination:
            %Command(box.port, '5 dim t');
            %Command(box.port, '10 while 1 do');
            % Print a character - Send 1 Byte to host computer:
            %Command(box.port, '15 let t = seconds % 10');
            %Command(box.port, '20 print t');
            % Sleep for 500 microseconds...
            %        Command(box.port, '30 sleep 1 ms');
            % Print a character - Send 1 Byte to host computer:
            %        Command(box.port, '35 let t = seconds % 10');
            %        Command(box.port, '40 print t');
            % Sleep for 500 microseconds...
            %        Command(box.port, '50 sleep 1 ms');
            % Then repeat.
            %Command(box.port, '60 endwhile');
        end
        
        % Print out whole program as debug output:
        Command(box.port, 'list');

        % Wait a second, then drain the input buffers:
        WaitSecs(1);
        while IOPort('BytesAvailable', box.port) > 0
            char(IOPort('Read', box.port))
        end
        
        % Box is programmed and ready. Start the program and thereby data
        % streaming.
        Command(box.port, 'run');
        
        % Let a second of data junk accumulate to simulate the real
        % situation where the box streams uncontrollably after power-up...
        WaitSecs(2);
        
        % Set input filter to discard CR and LF characters as well as
        % redundant data and attach a total streamcount tag and dT tag of
        % 32 bit size to each read datum, i.e., 2 * 32 bit = 8 bytes:
        IOPort('ConfigureSerialPort', box.port, 'ReadFilterFlags=3');
    else
        % Set input filter to discard redundant data and attach a total
        % streamcount tag and dT tag of 32 bit size to each read datum,
        % i.e., 2 * 32 bit = 8 bytes:
        IOPort('ConfigureSerialPort', box.port, 'ReadFilterFlags=1');
    end
    
    if box.type == 3
        % PST box: Disable streaming from box and switch off all lamps:
        % Bit 7 = 128 -> Enable/Disable streaming.
        % Bit 6 =  64 -> Lower bits control lamp state.
        % Bit 5 =  32 -> Enable/Disable button queries.
        % Bit 0-4 = 1-16 = Enable/Disable Lamp 0-4.
        IOPort('Write', box.port, uint8(64));
        
        % Wait extra 250 msecs for box to calm down:
        WaitSecs(0.25);
        
        % Now we can reasonably assume that no new data will arrive from
        % box and no data is in-flight. Whatever's there should be in the
        % OS receive buffers.
    end

    if box.type == 2
        % CMU box: Disable serial port receiver to stop the data stream
        % from box at the receiving side, as this braindead box can't be
        % instructed to stop streaming by itself. This call will only work
        % on OS/X and Linux and even there only on a small subset of serial
        % port hardware :-(
        IOPort('ConfigureSerialPort', box.port, 'ReceiverEnable=0');

        % Wait extra 250 msecs for box to calm down:
        WaitSecs(0.25);
        
        % Now we can at least pray that no new data will arrive from
        % box and no data is in-flight. Whatever's there should be in the
        % OS receive buffers.
    end

    % Preheat GetSecs:
    GetSecs;

    % Purge all input and output buffers:
    IOPort('Purge', box.port);
        
    % Now try to get rid of all of it:
    while IOPort('BytesAvailable', box.port) > 0
        IOPort('Read', box.port, 0);
    end
    
    if box.type == 3 || box.type == 2
        % Calibrate inter-byte-interval:
        if box.type == 3
            % PST box: Enable streaming:
            IOPort('Write', box.port, uint8(128+64+32+1+4+16), 0);
        else
            % CMU box: Enable receiver:
            IOPort('ConfigureSerialPort', box.port, 'ReceiverEnable=1');
        end

        % Wait until 8000 bytes have arrived, read them and timestamp:
        when1 = GetSecs;
        while IOPort('BytesAvailable', box.port) == 0
            % Still waiting for 1st byte...
            when1 = GetSecs;
        end
        
        nscanned = 0;
        while nscanned < 8000        
            % Wait for final byte: We must not request more than 1000 Bytes
            % per request, otherwise we'll get in trouble on OS/X!
            [dummy, when2] = IOPort('Read', box.port, 1, 1000);
            nscanned = nscanned + length(dummy);
        end
        
        if box.type == 3
            % PST: Stop streaming from box:
            IOPort('Write', box.port, uint8(64+2+8));
        else
            % CMU: Disable receiver:
            IOPort('ConfigureSerialPort', box.port, 'ReceiverEnable=0');            
        end
        
        % Wait extra 250 msecs for box to calm down:
        WaitSecs(0.25);

        % Purge all input and output buffers:
        IOPort('Purge', box.port);

        % Now try to get rid of all of it:
        while IOPort('BytesAvailable', box.port) > 0
            IOPort('Read', box.port, 0);
        end
        
        % All buffer drained, no more data. Compute real byteduration of a
        % single samplebyte:
        box.dt = (when2 - when1) / nscanned;

        if box.type == 2
            % If this is Windows, or if the results seem implausible for the
            % CMU box, ie., more than 0.1 msecs off the expected 1.0 msec, then
            % we simply hard-code dt to 1.0 msec.
            if IsWin || (box.dt < 0.0009) || (box.dt > 0.0011)
                % CMU Box: Assume 1 msec per sample:
                box.dt = 1/1000;
            end
        end
    else
        % Don't use interbyte interval on non CMU and non PST:
        box.dt = 0;
    end

    % Setup Async-Reads with blocking background read behaviour. This also
    % affects regular Sync-Reads by avoiding polling for the first byte!
    % This will cause the background async reader thread to block in the
    % serial port read function until data is available, then fetch,
    % process, timestamp and enqueue it, then repeat the cycle. No cpu
    % resources are wasted polling and the accuracy of the timestamps is
    % only limited by the scheduling latency of the operating system:
    IOPort('ConfigureSerialPort', box.port, 'BlockingBackgroundRead=1');

    if box.type == 3
        % PST box: Enable streaming from box:
        % Technically only bit 7 (128) controls streaming and its speed, the
        % rest controls lamp state (1 = Lamp dark, 0 = Lamp on) as
        % described above. A setting of 255 will enable streaming (128) and will
        % switch on all lamps due to 64 + 32 + 16 + 8 + 4 + 2 + 1:
        IOPort('Write', box.port, uint8(255), 0);
    end
    
    if box.type == 2
        % CMU box: Enable serial port receiver to start the data stream
        % from box at the receiving side, as this braindead box can't be
        % instructed to start streaming by itself. This call will only work
        % on OS/X and Linux and even there only on a small subset of serial
        % port hardware :-(
        IOPort('ConfigureSerialPort', box.port, 'ReceiverEnable=1');
    end
    
    % Fetch 1st sample synchronously, so we don't get a skewed box.baseTime
    % due to startup latencies of the async reader thread:
    [box.oldState, box.baseTime, box.olderr] = IOPort('Read', box.port, 1, 1);
    
    % Make sure we have a non-empty oldState:
    if isempty(box.oldState)
        box.oldState = -1;
    end
    
    % Start background read operation, try to fetch and timestamp data at a
    % granularity of 9 Byte -- Each single status byte from the box gets
    % timestamped individually, and a 4 Byte streamcount of read bytes gets
    % attached, as well as a 4 byte timedelta in microseconds, resulting in
    % 9 Bytes of data for each single non-redundant byte of data from the
    % box:
    IOPort('ConfigureSerialPort', box.port, 'StartBackgroundRead=9');
    
    % CAUTION: As soon as StartBackgroundRead has been called, we should
    % avoid calling IOPort('Write') or IOPort('Purge') etc. on MS-Windows.
    % Apparently access to serial ports is not thread-safe and we can get
    % into race-conditions if trying to access the port while our
    % background thread accesses as well. Only 'Read', 'BytesAvailable',
    % and a subset of 'ConfigureSerialPort' commands are safe on Windoze.
    
    % baseTime is first oldTime:
    box.oldTime = box.baseTime;
    
    % Reset trouble counter:
    box.tTrouble = 0;
    
    % Assign box struct to internal boxes struct array at slot box.port+1:
    boxes{box.port + 1} = box;
    
    % Return box handle:
    varargout{1} = box.port + 1;
    
    % Perform one dummy read to discard first event: It is always invalid.
    IOPort('Read', box.port, 1, 9);
    
    return;
end

% Close box:
if strcmpi(cmd, 'Close')
    if isempty(handle)
        error('CMUBox: Close: No "handle" for box to close provided!');
    end
    
    if length(handle) ~= 1
        error('CMUBox: Close: Passed argument is not a valid response box "handle"!');
    end
    
    if length(boxes) < handle
        error('CMUBox: Close: Passed argument is not a valid response box "handle"!');
    end

    % Retrieve box struct:
    box = boxes{handle};
    
    if isempty(box.useBitwhacker)
        error('CMUBox: Close: Passed argument is not a valid response box "handle"!');
    end

    %    IOPort('Verbosity', 10);
    
    % Stop async read, release queues:
    IOPort('ConfigureSerialPort', box.port, 'StopBackgroundRead');

    % Emulated box?
    if box.useBitwhacker
        % Stop streaming: Send a CTRL+C control character (ascii code 3).
        % This will abort the program running on the Bitwhacker:
        IOPort('Write', box.port, char(3));
        WaitSecs(0.2);

        % Send reset signal to StickOS:
        % Command(box.port, 'reset');
    end
    
    if box.type == 3
        % PST box: Disable streaming from box and all lamps:
        IOPort('Write', box.port, uint8(64));
        
        % Wait extra 250 msecs for box to calm down:
        WaitSecs(0.25);
        
        % For the fun of it, a little scan of the lights as a goodbye:
        for j=1:4
            i = 1;
            while i < 32
                IOPort('Write', box.port, uint8(64 + i));
                i = i * 2;
                WaitSecs(0.1);
            end

            i = 16;
            while i >= 1
                IOPort('Write', box.port, uint8(64 + i));
                i = i / 2;
                WaitSecs(0.1);
            end
        end
        
        % Flush lights a last time:
        IOPort('Write', box.port, uint8(64 + 32 + 16 + 8 + 4 + 2 + 1));
        WaitSecs(0.4);
        
        % Turn off lights:
        IOPort('Write', box.port, uint8(64));
        WaitSecs(0.4);

        % Everything off:
        IOPort('Write', box.port, uint8(0));
        WaitSecs(0.25);

        % Now we can reasonably assume that no new data will arrive from
        % box and no data is in-flight. Whatever is there should be in the
        % OS receive buffers.
    end
    
    % Close down driver and port:
    IOPort('Close', box.port);
    
    % Invalidate entry in boxes() cell array to make handle invalid:
    box.useBitwhacker = [];
    boxes{handle} = box;
    
    % Done.
    return;
end

if strcmpi(cmd, 'Status')
    if isempty(handle)
        error('CMUBox: Status: No "handle" for box to close provided!');
    end
    
    if length(handle) ~= 1
        error('CMUBox: Status: Passed argument is not a valid response box "handle"!');
    end
    
    if length(boxes) < handle
        error('CMUBox: Status: Passed argument is not a valid response box "handle"!');
    end

    % Retrieve box struct:
    box = boxes{handle};
    
    if isempty(box.useBitwhacker)
        error('CMUBox: Status: Passed argument is not a valid response box "handle"!');
    end
    
    % Return box status struct:
    varargout{1} = box;
    
    return;
end

% Invalid command!
error('CMUBox: Invalid or unknown command specified!');

% End of main function:
end

% Helper function for sending commands to Bitwhacker:
function Command(p, cmd)
    IOPort('Write', p, [cmd char(13)]);
    % WaitSecs(0.010);
end
