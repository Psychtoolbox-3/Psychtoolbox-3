function serialnumStr = PR670getserialnumber
% PR670getserialnumber - Gets the PR-670 serial number.
%
% Syntax:
% serialnumStr = PR670getserialnumber
%
% Description:
% Get the serial number from the PR-670.  Sometimes it's useful to be able
% to find out  which of the instruments in the lab was used for any
% particular measurement.
%
% Output:
% serialnumStr (1xN char) - The serial number as a char array.

global g_serialPort;

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Initialize
timeout = 30;

% Flush the buffers.
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR670read;
end

% Send the command to return the serial number.
PR670write('D110');

% Get at least one character
waited = 0;
inStr = [];
while isempty(inStr) && (waited < timeout)
	WaitSecs(1);
	waited = waited + 1;
	inStr = PR670read;
end
if waited == timeout
	error('No response from meter');
end

% Pick up entire buffer.
readStr = inStr;
while ~isempty(inStr)
	inStr = PR670read;
	readStr = [readStr inStr]; %#ok<AGROW>
end

% Parse return.  This may contain training blanks.
A = textscan(readStr, '%d,%s');
readStr = A{2}{1};
serialnumStr = readStr(1:8);
