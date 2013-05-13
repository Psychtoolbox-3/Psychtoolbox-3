function [xyz,qual] = PR655measxyz
% [xyz,qual] = PR655measspd
%
% Make a measurement of the spectrum.
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

% See if we can sync to the source
% and set sync mode appropriately.
syncFreq = PR655getsyncfreq;
if ~isempty(syncFreq) && syncFreq ~= 0
	PR655write('SS1');
else
    PR655write('SS0');
    disp('Warning: Could not sync to source.');
end

% Make measurement and get the string;
readStr = PR655rawxyz(timeout);

% Check returned data.
idx = findstr(readStr, ',');
qual = str2num(readStr(1:idx(1)-1));

% Check for other error conditions
if qual == -1 || qual == 10
    %disp('Low light level during measurement');
    %disp('Setting returned value to zero');
    xyz = zeros(3,1);
elseif qual == 18 || qual == 0
    units = str2num(readStr(idx(1)+1:idx(2)-1));
    if units ~= 0
        error('Units not returned as cd/m2');
    end
    xyz = str2num(readStr(idx(2)+1:end))';
elseif (qual ~= 0)
    error('Bad return code %g from meter',qual);
end

return
