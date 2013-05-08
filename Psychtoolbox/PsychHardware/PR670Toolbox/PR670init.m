function retval = PR670init(portNumber)
% PR670init - Initialize the serial port to talk to the PR-670.
%
% Syntax:
% retval = PR670init
% retval = PR670init(portNumber)
%
% Description:
% Initializes the serial port for talking to colorimeter. Returns whatever
% character is sent by colorimeter
%
% Input:
% portNumber (1xN char) - Name of the serial port to use.
%
% Output:
% retval (1xN char) - Character data returned from the PR-670.

global g_serialPort;

if nargin == 0
    % This seems to be the default name on OS X. We don't know about other
    % operating systems defaults:
    portNumber = FindSerialPort('usbmodem', 1);
end

% Only open if we haven't already.
if isempty(g_serialPort)
    % IOPort has above port settings 9600 baud, no parity, 8 data bits,
    % 1 stopbit, no handshake (aka FlowControl=none) already as
    % built-in defaults, so no need to pass them:
    oldverbo = IOPort('Verbosity', 2);
    g_serialPort = IOPort('OpenSerialPort', portNumber, 'Lenient DontFlushOnWrite=1');
    IOPort('Verbosity', oldverbo);
end

% Tell the PR-670 to exit remote mode.  We do this to ensure a response from
% the device when we go into remote mode.  I've found it only responds the
% first time it's asked to go into remote mode.  There is no return code
% for this command.
PR670write('Q', 0);
pause(0.5);

% Put in remote mode.
PR670write('PHOTO', 0);

% Get the response.  Timeout after 10 seconds.
StartTime = GetSecs;
retval = [];
while isempty(retval) && (GetSecs-StartTime) < 10
    retval = PR670read;
end
