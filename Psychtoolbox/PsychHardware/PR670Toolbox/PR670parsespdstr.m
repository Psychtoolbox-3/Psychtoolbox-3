function spd = PR670parsespdstr(readStr, S)
% PR670parsespdstr - Parses the  spectral power distribution string returned by the PR-670.
%
% Syntax:
% spd = PR670parsespdstr(readStr)
% spd = PR670parsespdstr(readStr, S);
%
% Description:
% Parse the spectral power distribution string returned by the PR-670.  The
% results are splined to the desired wavelength sampling defined by the 'S'
% input parameter.
%
% Input:
% readStr (1xN char) - Raw data returned from a meter measurment.
% S (1x3) - Wavelength sampling.  Default: [380 5 81]

if nargin < 2 || isempty(S)
	S = [380 5 81];
end

% Split up the data delimited by newline characters.
C = textscan(readStr, '%s', 'Delimiter', '\n');
C = C{1};

% Loop over all the data lines and pull out the spectral info.  Note that
% the first line doesn't contain wavelength data so we toss it.
C = C(2:end);
numWavelengths = length(C);
if (numWavelengths ~= 201)
    error('Unexpected native number of wavelength samples for PR-670');
end
spd = zeros(numWavelengths, 1);
wavelengths = spd;
for i = 1:length(C)
	% Parse the wavelength measurement line.  The first element of the
	% returned cell array will be the wavelength, the second element the
	% measurement.
	D = textscan(C{i}, '%d,%f');
	
	wavelengths(i) = D{1};
	spd(i) = D{2};
end

% Convert to our units standard.
spd = 2 * spd;
if (wavelengths(2)-wavelengths(1) ~= 2)
    error('Unexpected native wavelength sampling for PR-670');
end

% Spline to desired wavelength sampling.
startWavelength = wavelengths(1);
spacing = wavelengths(2) - wavelengths(1);
S0 = [startWavelength, spacing, numWavelengths];
spd = SplineSpd(S0, spd, S);
