function peaks = FindSpectralPeaks(spectrum,S,noiseLevel)
% peaks = FindSpectralPeaks(spectrum,[S],[noiseLevel])
%
% Find the local peaks in a spectrum.
% Uses the heuristic of requiring increasing
% and decreasing in neighborhood of peak to
% get rid of noise.
%
% Parameter noiseLevel (default 0) causes
% the routine to ignore peaks lower than
% noiseLevel*maxValue where maxValue is the
% largest measurement in the spectrum.
%
% 1/6/96		dhb		Wrote it.
% 5/17/99   dhb   Added  noiseLevel parameter.

if (nargin < 2 || isempty(S))
	S = [380 5 81];
end
wls = MakeItWls(S);
if (nargin < 3 || isempty(noiseLevel))
	noiseLevel = 0;
end

indices = [];
[nWls,nil] = size(spectrum);
noiseThreshold = noiseLevel*max(spectrum);
for i = 3:nWls-2
	if (spectrum(i-2) < spectrum(i-1) && ...
		  spectrum(i-1) < spectrum(i) && ...
			spectrum(i+1) < spectrum(i) && ...
			spectrum(i+2) < spectrum(i+1) && ...
			spectrum(i) > noiseThreshold)
		indices = [indices ; i];
	end
end


peaks = indices;
[nPeaks,nil] = size(indices);
for i = 1:nPeaks
	peaks(i) = wls(indices(i));
end

