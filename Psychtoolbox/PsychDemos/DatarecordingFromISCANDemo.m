function DatarecordingFromISCANDemo(portSpec, sampleFreq, maxComponents, baudRate, specialSettings)
% Template for asynchronous data collection and timestamping from serial port.
%
% DatarecordingFromSerialPortDemo([portSpec=auto-detect][, sampleFreq=240][maxComponents=2][, baudRate=115200][, specialSettings=None])
%
% This demo shows how to perform efficient data recording of data from an
% ISCAN eyetracker which is connected to the serial port or a USB-Serial
% converter. The device is assumed to stream data in packets of at most
% 'maxComponents' gaze sample components in ASCII streaming mode.
%
% The demo connects to the first found serial port, or optionally the port
% given by the 'portSpec' parameter, e.g., 'COM5'. It connects at a
% baudrate of 'baudRate' Baud, by default without flow-control, with 8
% databits, 1 stopbit and no parity, but you can set arbitrary settings via
% the optional 'specialSettings' string (see IOPort OpenSerialPort?
% online-help for possible parameters).
%
% Then it allocates receivebufferspace for up to 1 hour of uninterrupted
% recording, then starts background recording of data.
%
% Datapackets can be read out in realtime, as demonstrated in the main
% while loop here, or offline at end of a session. Each packet is padded to
% be 'maxReadQuantum' bytes in size (zeros are added for shorter packets).
% Each packet comes with a GetSecs() timestamp of when the first byte of a
% packet was received.
%
% Of course you'll need to understand the code of this demo and then
% customize it for your needs. We don't know what kind of weird stuff
% you're going to connect, but the default parameter settings may work on
% some simple devices like eyetrackers.
%
% For accurate timestamping and data reception with low latency, make sure
% that you've configured your serial ports properly. Search the
% Psychtoolbox forum for posts on that topic, e.g., message 9873.
%
% The following setup has to be performed in the ISCAN DQW software,
% citing ISCAN software documentation:
%
% To enable the serial output capability in DQW: In the DQW "System
% Configuration" dialog, select the "Serial I/O" tab. Of the available COM
% channels, specify the COM port to be used as "Data I/O", the baud rate to
% match the receiving system, and select "Output" as active. Click on "OK"
% to register changes to the system configuration. Click on the "Options"
% button until the "Auxiliary Output Controls" panel appears as the lower
% left DQW screen. Then click on the "Serial" tab if this is not already in
% the foreground. Choose "Raw Binary" or "Raw ASCII" as the output type
% depending on your preference. Two parameter banks with six parameters
% each, for a total of up to 12 parameters may be output with each data
% sample. To begin, select "Param Bank" 1 and fill the channels (01 -> 06)
% with the desired parameters by selecting from the pop up list activated
% by clicking on the "..." button to the right of each channel. Fill the
% channels in order and be sure that unused channels have parameter
% "................" selected to minimize data transmission time. If more
% than 6 parameters are desired, select "Param Bank" 2 and proceed as above
% to fill in channels 07 -> 12. With the above steps completed, whenever
% "Track Active" is enabled a data packet will be output via the specified
% COM port for each system data sample. The data packet will consist of
% each of the non-blank channel data values for the parameters specified in
% channel order (01 -> 12).
%
% To enable the serial input capability in DQW: In the DQW "System
% Configuration" dialog, select the "Serial I/O" tab. Of the available COM
% channels, specify the COM port to be used as "Data I/O", the baud rate to
% match the remote data sending system, and select "Input" as active. Click
% on "OK" to register changes to the system configuration. With the above
% steps completed, the DQW system is ready to receive serial commands from
% a remote source. All DQW remote serial commands are single bytes between
% 128 and 255 (80 -> FF Hex). Whenever "Track Active" is enabled, the DQW
% system can receive and record any single serial byte between 0 and 127
% (00 and 7F Hex) as a synchronizing data marker. These bytes are
% represented as parameter "SerIn0" wherever parameter selections are made
% for display or recording in the DQW program. To trigger data recording
% via the serial port, in the DQW "System Configuration" dialog, select the
% ?Recording? tab and set the recording ?Trigger? to ?Serial.? 
%
% Thanks to Sebastian Moeller for finding this out.
%

% History:
% 16.9.09  mk  Written, derived from DataRecordingFromSerialPortDemo.

% joker variable: Ignore this! Only for MK's internal testing:
joker = '';

% MK: Only for testing at MK's site: joker = 'Lenient';

% maxReadQuantum provided?
if ~exist('maxComponents', 'var')
    maxComponents = [];
end

if isempty(maxComponents)
    % Assume that a single datapacket is no more than 16 Bytes, including
    % any terminator bytes:
    maxComponents = 2;
end

% Each component is 7 Bytes (space + xxx.x digits + tab delimiter), and
% there is a CR+LF combo at the end with the lineTerminator being the last
% entry:
maxReadQuantum = maxComponents * 7 + 2;

% Line terminator for ISCAN is ASCII code 10, aka LF, aka NL aka newline/linefeed:
lineTerminator = 10;

% sampleFreq provided?
if ~exist('sampleFreq', 'var')
    sampleFreq = [];
end

if isempty(sampleFreq)
    % Choose an expected sampling frequency for incoming data packets of
    % 240 Hz for the ISCAN:
    sampleFreq = 240;
end

% Baudrate provided?
if ~exist('baudRate', 'var')
    baudRate = [];
end

if isempty(baudRate)
    % Choose an optimitic default of 115.2 KBaud:
    baudRate = 115200;
end

% Any serial port specification provided?
if ~exist('portSpec', 'var')
    portSpec = [];
end

if isempty(portSpec)
    % Try to auto-detect a suitable serial port:
    portSpec = FindSerialPort([], 1);
end

% Any serial specialSettings provided?
if ~exist('specialSettings', 'var')
    % Set to empty aka "use defaults" which is 8-N-1 encoding, no flow
    % control...
    specialSettings = [];
end

% Compute maximum input buffer size for 1 hour worth of datapackets coming
% in at a expected sampleFreq Hz with a size of at most maxReadQuantum Bytes
% each:
InputBufferSize = maxReadQuantum * sampleFreq * 3600;

% Assign an interbyte readtimeout which is either 15 seconds, or 10 times
% the expected time between consecutive datapackets at the given sampleFreq
% sampling frequency, whatever's higher. Could go higher or lower than
% this, but this seems a reasonable starter: Will give code and devices
% time to start streaming, but will prevent script from hanging longer than
% 15 seconds if something goes wrong with the connection:
readTimeout = max(10 * 1/sampleFreq, 15);

% HACK: Restrict maximum timeout to 21 seconds. This is needed on Macintosh
% computers, because at least OS/X 10.4.11 seems to have a bug which can
% cause the driver to hang when trying to stop at the end of a session if
% the timeout value is set higher than 21 seconds!
readTimeout = min(readTimeout, 21);

% Assemble initial configuration string for opening the port with
% reasonable settings: Given special settings, baudrate, inputbuffersize.
% Also set the special delimiter character code 'lineTerminator' that
% signals the end of a valid data packet:
portSettings = sprintf('%s %s BaudRate=%i InputBufferSize=%i Terminator=%i ReceiveTimeout=%f', joker, specialSettings, baudRate, InputBufferSize, lineTerminator, readTimeout );

% Open port portSpec with portSettings, return handle:
myport = IOPort('OpenSerialPort', portSpec, portSettings);

fprintf('Link online: Hit a key on keyboard to start data recording, after that hit any key to finish data collection.\n');
KbStrokeWait;

% ---- Here you'd put any IOPort setup calls, e.g., write and read commands
% to setup your device, enable streaming of data etc...

iscan_start_code = char(hex2dec('84'));
iscan_stop_code = char(hex2dec('88'));

% Stop tracker, if not already stopped:
IOPort('Write', myport, iscan_stop_code);

% Wait a bit...
WaitSecs(1);

% Purge data buffers, read and discard remaining junk:
IOPort('Purge', myport);

WaitSecs(0.5);

IOPort('Read', myport);

% Really empty? Hope so.
if IOPort('BytesAvailable', myport) > 0
    warning('Trouble ahead!! Check your eyetracker!'); %#ok<WNTAG>
end

% All buffers empty, device idle, ready to go:

% Configure background read.
asyncSetup = sprintf('%s BlockingBackgroundRead=1 ReadFilterFlags=4', joker);
IOPort('ConfigureSerialPort', myport, asyncSetup);

% Start tracker:
IOPort('Write', myport, iscan_start_code);

% ---- End of device specific setup ----

% Start asynchronous background data collection and timestamping. Use
% blocking mode for reading data -- easier on the system:
asyncSetup = sprintf('%s StartBackgroundRead=%i', joker, maxReadQuantum);
IOPort('ConfigureSerialPort', myport, asyncSetup);

% Data collection started: From now on, the driver will read data from the
% serial port, byte by byte. Whenever either maxReadQuantum bytes have
% accumulated or the packet delimiter code 'lineTerminator' is detected,
% the driver assumes a packet is complete. If less than 'maxReadQuantum'
% bytes have been received before the 'lineTerminator' is detected, the
% driver will pad the remaining bytes with zeros. --> Any packet will
% always be exactly 'maxReadQuantum' bytes in size, zero-padded or not.
%
% Reception of the first byte of a new packet is timestamped in GetSecs()
% time, and the later IOPort('Read') calls will return those timestamps...

% ---- From here on, only a limited set of IOPOrt commands is allowed if
% you are working on MS-Windows. On OS/X and Linux no such restrictions
% apply afaik.


% Alloc timestamp array for 1 hour worth of samples at 'sampleFreq' Hz:
tpkt = zeros(1, sampleFreq * 3600);
count = 0;

% Fake experiment loop: Run until any key is pressed:
while ~KbCheck
    % Wait blocking for a new data packet of 'maxReadQuantum' bytes from
    % the serial port, then return the packet data as uint8's plus the
    % GetSecs receive timestamp 'treceived' of the start of each packet:
    [pktdata, treceived] = IOPort('Read', myport, 1, maxReadQuantum);
    
    % Empty data packet?
    if ~any(pktdata)
        % All zeroes: An empty packet due to error? Skip it.
        fprintf('Empty packet received (all zeros) --> SKIPPED\n');
        % Back to start of while loop:
        continue;
    end
    
    % Some data available to process. For the sake of demonstration assume
    % the pktdata actually encodes an ASCII text string. Convert it to a
    % string, print it and its timestamp. We deblank() to strip trailing
    % zero bytes:
    pktstring = deblank(char(pktdata));
    fprintf('ReceiveTime: %f seconds GetSecs time: Datastring = %s\n', treceived, pktstring);
    
    % Store timestamp for the fun of it:
    count = count + 1;
    tpkt(count) = treceived;
    
    % HERE YOU'D DO something useful with the data...
    
    
    % Next loop iteration...
end

% We'd like to stop data collection here. However there may be still data
% queued up in the receive buffers, so we will now try to get all data that
% has been received up to this point in time:
tEnd = GetSecs;

fprintf('TRIAL LOOP STOPPED AT t = %f seconds. Now fetching pending data up to that point...\n', tEnd);

% Fetch all pending data that has been received up to systemtime tEnd:
while treceived < tEnd
    % Same as above, but now a non-blocking read (flag == 0), ie., once the
    % end of the dataqueue is reached, simply return an empty 'pktdata'
    % variable immediately:
    [pktdata, treceived] = IOPort('Read', myport, 0, maxReadQuantum);
    
    % Empty data packet?
    if isempty(pktdata) || ~any(pktdata)
        % All zeroes or an empty packet. We stop here.
        fprintf('Empty packet received (all zeros) --> END OF DATA FETCH\n');
        % Break out of loop:
        break;
    end
    
    % Some data available to process. For the sake of demonstration assume
    % the pktdata actually encodes an ASCII text string. Convert it to a
    % string, print it and its timestamp. We deblank() to strip trailing
    % zero bytes:
    pktstring = deblank(char(pktdata));
    fprintf('ReceiveTime: %f seconds GetSecs time: Datastring = %s\n', treceived, pktstring);
    
    % Store timestamp for the fun of it:
    count = count + 1;
    tpkt(count) = treceived;

    % Next iteration...
end
    
% All data up to time tEnd fetched. There may be more data lingering in the
% queue but we don't care...


% End of data recording. Stop background read operation. This will discard
% any remaining data in the receive buffers. Stopping is a bit tricky if
% the source device doesn't transmit anymore. Then this may hang until the
% ReceiveTimeout timeout value expires, ie., up to 'readTimeout' seconds.
%
% Therefore it is better to only stop data transmission from the sending
% device after we have executed the following stop commands.
fprintf('STOPPING DATA COLLECTION: This may take up to %f seconds...\n', readTimeout);
% This line may or may not help if you experience hangs at this place: IOPort('ConfigureSerialPort', myport, [joker ' BlockingBackgroundRead=0']);
IOPort('ConfigureSerialPort', myport, [joker ' StopBackgroundRead']);

% Now the driver has discarded the data and is in synchronous manual mode
% of operation again. You could add any kind of IOPort commands to shut
% down your serial device now...
fprintf('SHUTTING DOWN DEVICE\n');

% Stop tracker:
IOPort('Write', myport, iscan_stop_code);

% Wait a bit...
WaitSecs(1);

% Close port and driver:
fprintf('CLOSING SERIAL PORT\n');
IOPort('Close', myport);
fprintf('DONE.\n');

% Plot delta between consecutive received packets for the fun of it:
tpkt = tpkt(1:count);

% Deltas in msecs:
tpkt = 1000 * diff(tpkt);

% Print mean delay and sample frequency:
fprintf('Mean delay between packets is %f msecs, i.e., %f Hz real datarate vs. expected %f Hz.\n', mean(tpkt), 1000/mean(tpkt), sampleFreq);

% Plot it:
close all;
plot(tpkt);

fprintf('Done. Bye!\n');

return;
