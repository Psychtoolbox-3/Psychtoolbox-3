function cornealIrradiance_PowerPerArea = RadianceAndDegrees2ToCornIrradiance(radiance_PowerPerSrArea,stimulusAreaDegrees2)
% cornealIrradiance_PowerPerArea = RadianceAndDegrees2ToCornIrradiance(radiance_PowerPerSrArea,stimulusAreaDegrees2)
%
% Convert the radiance of a stimulus to corneal irradiance, given that we know the area of the stimulus in degrees2.
% The routine assumes that the stimulus is rectangular with linear subtense sqrt(stimulusAreaDegrees2).
%
% Light power can be in your favorite units (Watts, quanta/sec) as can distance (m, cm, mm).  The units for
% area in the returned irradiance match those used for area in the passed radiance.
%
% So, if radiance is in Watts/[cm2-sr] then distance needs to be in cm and irradiance will be in Watts/cm2.
%
% The derivation also assumes the small angle approximation simulusSizeUnits = stimulusSizeRadians*stimulusDistanceUnits,
% where units are the relavant units of length.  Although we don't have stimulusSizeUnits and stimulusDistanceUnits,
% these turn out to cancel out under the small angle approximation.
%
% See also: RadianceAndDistanceAreaToCornIrradiance
%
% 2/20/13  dhb  Wrote it.

% Convert area in degrees squared to linear angluar subtense in radians.
stimulusSizeDegrees = sqrt(stimulusAreaDegrees2);
stimulusSizeRadians = deg2rad(stimulusSizeDegrees);

% Note that the following is true in the small angle approximation.
%   stimulusSize/stimulusDistance = stimulusSizeRadians
% where size and distance are in units like mm, cm, m etc.
%
% Note also that 
%   cornealIrradiance_PowerPerArea = radiance_PowerPerSrArea*((stimulusSize/stimulusDistance)^2)
% See comments in RadianceAndDistanceAreaToCornIrradiance for the derivation of this.
%
% These two observations yield the answer
cornealIrradiance_PowerPerArea = radiance_PowerPerSrArea*(stimulusSizeRadians^2);

