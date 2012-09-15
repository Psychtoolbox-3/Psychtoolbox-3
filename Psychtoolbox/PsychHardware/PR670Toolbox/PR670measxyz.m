function [xyz, qual] = PR670measxyz
% PR670measxyz - Make an XYZ measurement with the PR-670.
%
% Syntax:
% [xyz, qual] = PR670measxyz
%
% Output:
% xyz (3x1) - Measured XYZ values.
% qual (scalar) - Measurement quality code.  Will be 0 for successful
%     measurements, -8 if too little light.

global g_serialPort;

% Check for initialization
if isempty(g_serialPort)
   error('Meter has not been initialized.');
end

% Set the response timeout in seconds.
timeout = 30;

% See if we can sync to the source and set sync mode appropriately.
syncFreq = PR670getsyncfreq;
if ~isempty(syncFreq) && syncFreq ~= 0
	PR670write('SS1');
else
	disp('Warning: Could not sync to source.');
    PR670write('SS0');
end

% Make measurement and get the string;
readStr = PR670rawxyz(timeout);

% Extract the returned data.
C = textscan(readStr, '%d,%d,%f,%f,%f');

% Check returned data.
qual = C{1};

switch qual
	% Measurement OK
	case 0
		% Make sure the units look OK.
		units = C{2};
		assert(units == 0, 'Units not returned as cd/m2');
		
		% Extract the XYZ values.
		X = C{3};
		Y = C{4};
		Z = C{5};
		
		xyz = [X;Y;Z];
		
	% Low light
	case -8
		xyz = zeros(3,1);
		
	otherwise
		error('Bad return code %g from meter', qual);
end
