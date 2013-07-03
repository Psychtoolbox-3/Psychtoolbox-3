function radiance_PowerPerSrArea = CornIrradianceAndDegrees2ToRadiance(cornealIrradiance_PowerPerArea,stimulusAreaDegrees2)
% radiance_PowerPerSrArea = CornIrradianceAndDegrees2ToRadiance(cornealIrradiance_PowerPerArea,stimulusAreaDegrees2)
%
% Convert the corneal irradiance of a stimulus to radiance, given that we know the area of the stimulus in degrees2.
% The routine assumes that the stimulus is rectangular with linear subtense sqrt(stimulusAreaDegrees2).
%
% Light power can be in your favorite units (Watts, quanta/sec) as can distance (m, cm, mm).  The units for
% area in the returned radiance match those used for area in the passed irradiance.  So, if irrradiance is in Watts/cm2
% the radiance will be in Watts/[cm2-sr].
%
% The derivation assumes the small angle approximation simulusSizeUnits = stimulusSizeRadians*stimulusDistanceUnits,
% where units are the relavant units of length.  Although we don't have stimulusSizeUnits and stimulusDistanceUnits,
% these turn out to cancel out under the small angle approximation.
%
% See also: RadianceAndDistanceAreaToCornIrradiance, RadianceAndDegrees2ToCornIrradiance
%
% 2/22/13  dhb  Wrote it.

% Convert area in degrees squared to linear angluar subtense in radians.
stimulusSizeDegrees = sqrt(stimulusAreaDegrees2);
stimulusSizeRadians = deg2rad(stimulusSizeDegrees);

% This routine just inverts the simple relation derived in the comments to
%   RadianceAndDegrees2ToCornIrradiance
radiance_PowerPerSrArea = cornealIrradiance_PowerPerArea/(stimulusSizeRadians^2);

