function PR670setsyncfreq(syncFreq)
% PR670setsyncfreq - Sets the sync frequency for the light source.
%
% Syntax:
% PR670setsyncfreq(syncFreq)
%
% Input:
% syncFreq (scalar) - The sync frequency for source.  Can be 1 of 3 options:
%     1. 0 - No sync
%     2. 1 - Auto sync
%     3. 20-400 - Frequency in Hz
%
% Notes:
% See also PR670getsyncfreq

global g_serialPort

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Initialize
timeout = 30;

% Flushing buffers.
% fprintf('Flush\n');
dumpStr = '0';
while ~isempty(dumpStr)
	dumpStr = PR670read;
end

% Set the sync frequency.
switch syncFreq
	% Turn off sync mode.
	case 0
		PR670write('SS0');
		
		% Check the response.
		checkResponse(timeout, 'No response after SS0 command');
	
	% Auto sync.
	case 1
		PR670write('SS1');
		
		% Check the response.
		checkResponse(timeout, 'No response after SS1 command');
		
	% User defined frequency.
	otherwise
		% Force the frequency to an integer.  Not 100% sure if this is
		% required, but the documentation makes it look like only integers
		% are allowed.
		syncFreq = round(syncFreq);
		
		% Make sure the frequency is legit.
		assert(syncFreq >= 20 && syncFreq <= 400, 'PR670setsyncfreq:InvalidFreq', ...
			'Invalid sync frequency of %f, must be in the range [20,400]', ...
			syncFreq);
		
		% Tell the device we're specifying a user defined frequency.
		PR670write('SS3');
		
		% Check the response.
		checkResponse(timeout, 'No response after SS3 command');
		
		% Send the frequency over.
		PR670write(sprintf('SK%.3d', syncFreq));
		
		% Check the response.
		checkResponse(timeout, 'No response after SK command');
end


function checkResponse(timeout, timeoutString)
% checkResponse - Check for a PR-670 response and makes sure it's OK.

waited = 0;
inStr =[];
while isempty(inStr) && (waited < timeout)
	WaitSecs(1);
	waited = waited+1;
	inStr = PR670read;
end
if waited == timeout
	error(timeoutString);
end

% Pick up entire buffer.  This is the loop referred to above.
response = inStr;
while ~isempty(inStr)
	inStr = PR670read;
	response = [response inStr]; %#ok<AGROW>
end

% Parse the return and make sure we got a 0.  Any other value means an
% error occured.
qual = sscanf(response, '%d', 1);
if qual ~= 0
	error('PR670setsyncfreq:Received error code: %d', qual);
end
