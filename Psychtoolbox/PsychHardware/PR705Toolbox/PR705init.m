function response = PR705init(portString)
% PR705init - Initialize the serial port to talk to the PR-705.
%
% Syntax:
% retval = PR705init
% retval = PR705init(portString)
%
% Description:
% Initializes the serial port and establishes remote mode for the PR-705.
%
% Input:
% portString (1xN char) - Name of the serial port to use.
%
% Output:
% response (1xN char) - Character data returned from the PR-705.
%
% 11/29/12    zlb   Wrote it based on the PR670Toolbox. 

global g_serialPort g_useIOPort

if nargin == 0
    % Let FindSerialPort do its job
    portString = FindSerialPort('', g_useIOPort);
end

% Only open if we haven't already.
if isempty(g_serialPort)
    % IOPort has above port settings 9600 baud, no parity, 8 data bits,
    % 1 stopbit, no handshake (aka FlowControl=none) already as
    % built-in defaults, so no need to pass them:
    oldverbo = IOPort('Verbosity', 2);
    g_serialPort = IOPort('OpenSerialPort', portString, 'Lenient DontFlushOnWrite=1');
    IOPort('Verbosity', oldverbo);
end

PR705write('PR705');
response = PR705read(1, 16); % get the ' REMOTE MODE' response
