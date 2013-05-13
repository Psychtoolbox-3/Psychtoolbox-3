function retval = PR650init(portNumber, enableHandshaking)
% retval = PR650init(portNumber, [enableHandshaking])
% 
% Initialize serial port for talking to colorimeter.
% Returns whatever character is sent by colorimeter
%
% 'enableHandshaking' allows you to enable handshaking.  By default,
% handshaking is disabled.  To enable handshaking, set this value to 1 or
% true.
%
% 11/26/07    mpr   added timeout if nothing is returned within 10 seconds.
%
% In my experience, calling this function directly leads to poor performance
% (usually no communication is ever established).  You should find the function
% CMCheckInit in the PsychHardware folder, one folder up the tree from this one
% (which is presumably PR650Toolbox).  Calling this function from CMCheckInit
% should provide more reliable establishment of contact and hints on what to 
% try if contact fails.  -- MPR
 
global g_serialPort g_useIOPort;

if nargin == 1
    enableHandshaking = 0;
end

if ~g_useIOPort
    if enableHandshaking
        handshakeCode = 'h';
    else
        handshakeCode = 'n';
    end

    % Only open if we haven't already.
    if isempty(g_serialPort)
        SerialComm('open', portNumber, '9600,n,8,1');
        SerialComm('hshake', portNumber, handshakeCode);
        SerialComm('close', portNumber);
        WaitSecs(0.5);
        SerialComm('open', portNumber, '9600,n,8,1');
        SerialComm('hshake', portNumber, handshakeCode);
        g_serialPort = portNumber;
    end

    StartTime = GetSecs;
    % Send set backlight command to high level to check
    % whether we are talking to the meter.
    SerialComm('write', g_serialPort, ['b3' char(10)]);
    retval = [];
    while isempty(retval) && GetSecs-StartTime < 10
        retval = PR650serialread;
    end
end

if g_useIOPort
	if enableHandshaking
		handshakeCode = 'Lenient DontFlushOnWrite=1 FlowControl=Hardware ';
	else
		handshakeCode = 'Lenient DontFlushOnWrite=1 FlowControl=None ';
	end

	% Only open if we haven't already.
	if isempty(g_serialPort)
		oldverbo = IOPort('Verbosity', 2);
		hPort = IOPort('OpenSerialPort', portNumber, handshakeCode);
		IOPort('Close', hPort);
		WaitSecs(0.5);
		hPort = IOPort('OpenSerialPort', portNumber, handshakeCode);
		IOPort('Verbosity', oldverbo);
		g_serialPort = hPort;
	end

	StartTime = GetSecs;

	% Send set backlight command to high level to check
	% whether we are talking to the meter.
	IOPort('write', g_serialPort, ['b3' char(10)]);

	% Make sure the meter responds.
	retval = [];
	while isempty(retval) && GetSecs-StartTime < 10
		retval = PR650serialread;
	end
end
