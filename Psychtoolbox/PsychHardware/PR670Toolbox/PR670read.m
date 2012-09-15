function serialData = PR670read
% PR670read - Read data from the PR-670.
%
% Syntax:
% serialData = PR670read
%
% Description:
% Reads data chars from PR-670 until there is nothing left.  Returns an
% empty matrix if there is nothing to read.
% 
% Output:
% serialData (1xN char) - Data read from the PR-50.

global g_serialPort;

if isempty(g_serialPort)
	error('Meter has not been initialized.');
end

% Look for any data on the serial port.
serialData = char(IOPort('Read', g_serialPort));

% If data exists keep reading off the port until there's nothing left.
if ~isempty(serialData)
    tmpData = 1;
    while ~isempty(tmpData)
        WaitSecs(0.050);
        tmpData = char(IOPort('Read', g_serialPort));
        serialData = [serialData, tmpData]; %#ok<AGROW>
    end
end
