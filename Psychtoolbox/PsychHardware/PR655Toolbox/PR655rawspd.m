function readStr = PR655rawspd(timeout)
% readStr = PR655rawspd(timeout)
%
% Measure spd and return string.
% 
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
%

global g_serialPort;

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Flushing buffers.
% fprintf('Flush\n');
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR655read;
end

% Make measurement
% fprintf('Measure\n');
PR655write('M5')
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
