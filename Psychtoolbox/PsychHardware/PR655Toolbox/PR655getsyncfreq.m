function syncFreq = PR655getsyncfreq
% syncFreq = PR655getsyncfreq
%
% Measure sync frequency for source.  Returns
% empty if can't sync.
% 
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
%

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
	dumpStr = PR655read;
end


% Make measurement
PR655write('F')

StartTime = GetSecs;
waited = GetSecs-StartTime;
inStr = [];
while isempty(inStr) && waited < timeout
	inStr = PR655read;
    waited = GetSecs-StartTime;
end
if waited >= timeout
	error('No response after measure command');
end

readStr = inStr;

% Parse return
qual = -1;
[raw, count] = sscanf(readStr,'%f,%f',2);
if count == 2
	qual = raw(1);
	syncFreq = raw(2);
end

if qual ~= 0
	syncFreq = [];
end

return
