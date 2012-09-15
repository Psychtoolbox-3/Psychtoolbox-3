function readStr = PR670rawxyz(timeout)
% PR670rawxyz - Makes a raw XYZ measurement using the PR-670.
%
% Syntax:
% readStr = PR670rawxyz(timeout)
%
% Input:
% timeout (scalar) - The timeout period of the function in seconds.
%
% Output:
% readStr (1xN char) - The raw character output from the PR-670.

% Flush the buffers.
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR670read;
end

% Make measurement
PR655write('M2');
StartTime = GetSecs;
waited = GetSecs-StartTime;
inStr = [];
while isempty(inStr) && waited < timeout
	inStr = PR670read;
    waited = GetSecs-StartTime;
end
 
if waited == timeout
   error('Unable to get reading from radiometer');
else
	readStr = inStr;
end
