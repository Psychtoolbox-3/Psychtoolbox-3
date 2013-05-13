function readStr = PR655rawxyz(timeout)
% readStr = PR655rawxyz(timeout)
%
% Make a raw XYZ measurement.
% 
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
%

global g_serialPort;

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Flushing buffers.
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR655read;
end

% Make measurement
PR655write('M2')
StartTime = GetSecs;
waited = GetSecs-StartTime;
inStr = [];
while isempty(inStr) && waited < timeout
	inStr = PR655read;
    waited = GetSecs-StartTime;
end
 
if waited == timeout
   error('Unable to get reading from radiometer');
else
	readStr = inStr;
end

return
