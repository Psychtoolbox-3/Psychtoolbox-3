function [spd, qual] = PR670measspd(S, syncMode)
% PR670measspd - Make a measurement measurement.
%
% Syntax:
% [spd, qual] = PR670measspd(S)
% [spd, qual] = PR670measspd(S, syncMode)
%
% Input:
% S (1x3) - Desired wavelength sampling. Default: [380 5 81]
% syncMode (string) - Toggles syncronization with the light source.  Set
%     using the strings 'on' or 'off'.  Default: 'on'
%
% Output:
% spd (Mx1) - Light source spectrum.
% qual (scalar) - Measurement quality code.  See the PR-670 manual.

% Handle defaults
if nargin < 2 || isempty(syncMode)
    syncMode = 'on';
end

% Set wavelength sampling if passed.
if nargin < 1 || isempty(S)
   S = [380 5 81];
end

% Initialize
timeout = 300;

% See if we can sync to the source and set sync mode appropriately.
if strcmp(syncMode, 'on')
    syncFreq = PR670getsyncfreq;
    if ~isempty(syncFreq)
        PR670setsyncfreq(1);
    else
        PR670setsyncfreq(0);
    end
else
    PR670setsyncfreq(0);
end

% Get the meter response.
readStr = PR670rawspd(timeout);

% Extract the result code.
qual = sscanf(readStr, '%f', 1);
	 
switch qual
	% Measurement OK
	case 0
		spd = PR670parsespdstr(readStr, S);
		
	% Too dark
	case -8
		spd = zeros(S(3), 1);
		
	% Light source sync failure
	case {-1, -10}
		disp('Could not sync to source, turning off sync mode and remeasuring');
		PR670write('SS0');
		readStr = PR670rawspd(timeout);
		qual = sscanf(readStr,'%f',1);
		if qual == -8
			spd = zeros(S(3), 1);
		elseif qual == 0
			spd = PR670parsespdstr(readStr, S);
		else
			error('Received unhandled error code %d\n', qual);
		end
	
	otherwise
		error('Bad return code %g from meter', qual);
end
