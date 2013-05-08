function retval = PR655init(portNumber)
% retval = PR655init(portNumber, [enableHandshaking])
%
% Initialize serial port for talking to colorimeter.
% Returns whatever character is sent by colorimeter
%
% Per PhotoResearch, handshaking is disabled.
%
% 11/26/07    mpr   added timeout if nothing is returned within 10 seconds.
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
%
% It seems the iterative CMCheckInit method of initialization is not
% necessary with the PR655, so one run of this seems to do the trick.
% "usbmodem*" seems to cover every instance I've run into. Not sure about
% the prevalance of using usbmodem as a generic identifier, but if you can
% afford the PR655, you're probably on some more respectable form of internet
% connection. -TBC
%

global g_serialPort;

if nargin == 0
    % This seems to be the default name on OS/X. We don't know about other
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

% Put in Remote Mode --No [CR] after 'PHOTO'
rm = 'PHOTO';
for i = 1:5
    IOPort('write', g_serialPort, rm(i));
end

StartTime = GetSecs;
retval = [];
while isempty(retval) && GetSecs-StartTime < 10
    retval = PR655read;
end

return;
