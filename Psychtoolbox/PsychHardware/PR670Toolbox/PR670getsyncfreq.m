function [syncFreq, errorCode] = PR670getsyncfreq
% PR670getsyncfreq - Measures the frequency of the light source.
%
% Syntax:
% syncFreq = PR655getsyncfreq
%
% Output:
% syncFreq (scalar) - The frequency of the light source in Hz.  Is empty if
%     PR-670 encounters an error.
% errorCode (scalar) - The error code if an error occurs.  Is 0 if no
%     error.

global g_serialPort;

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Initialize
timeout = 30;

% Flushing buffers.
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR670read;
end

% Take a measurement.
PR670write('F');

% Get the response.
StartTime = GetSecs;
waited = GetSecs - StartTime;
inStr = [];
while isempty(inStr) && waited < timeout
	inStr = PR670read;
    waited = GetSecs - StartTime;
end
if waited >= timeout
	error('No response after measure command');
end

readStr = inStr;

% Parse the return string.
errorCode = -1;
[raw, count] = sscanf(readStr,'%f,%f',2);
switch count
	% Error occured
	case 1
		errorCode = raw(1);
		syncFreq = [];
	
	% No error
	case 2
		errorCode = raw(1);
		syncFreq = raw(2);
end
