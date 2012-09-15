function readStr = PR670rawspd(timeout)
% PR670rawspd - Takes an spd measurement and returns the results.
%
% Syntax:
% readStr = PR655rawspd(timeout)
%
% Input:
% timeout (scalar) - Timeout period in seconds.
%
% Output:
% readStr (1xN char) - The raw character array resulting from a
%     measurement.

% Flush the buffers.
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR670read;
end

% Make measurement
PR670write('M5')
StartTime = GetSecs;
waited = GetSecs - StartTime;
inStr = [];
while isempty(inStr) && waited < timeout
	inStr = PR670read;
    waited = GetSecs - StartTime;
end

if waited == timeout
   error('Unable to get reading from radiometer');
else
    readStr = inStr;
end
