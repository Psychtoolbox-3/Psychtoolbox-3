function [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousRetIrradiancePCWeightedValue(...
    S,radiance_WattsPerSrM2,weightingA,stimulusDurationSecs,eyeLengthMm)
% function [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousRetIrradiancePCWeightedValue(...
%   S,radiance_WattsPerSrM2,weightingA,stimulusDurationSecs,[eyeLengthMm])
%
% Compute the weighted aphakic (photochemical) retinal irradiance for Type 1 instruments as given on page 8, Table 2, 
% 5.4.1.3.a.
%
% Input spectrum is radiance in units of Watts/[sr-m2-wlinterval].
%
% Also return the exposure limit for this quantity.
%
% See page 6 for a definition of a Type 1 instrument.  As far as I can tell, the key
% criterion is that it doesn't put out more light that exceeds the Type 1 limits.
% 
% If the exposure time is longer than 2 hours the specified limits should be reduced by
% 1/exposureDuration in hours.  This routine implements that adjustment for its returned
% limit value.  It does not implement a further reduction of of the limit (by a factor of 2)
% specifed for microscopes and endoilluminators.
%
% The standard specifies that the passed radiance should be the highest averaged over
% an aperture of a specified size, where the size depends on the instrument.  This
% routine does not worry about that aspect.  The most conservative thing to do is
% to pass the highest localized power that will be presented.
%
% The standard specifies a pupil diameter (7 mm) to use for the conversion from radiance
% to retinal irradiance, but not an eye length.  We use 17 mm here by default.  You
% can override this by passing a different length in mm.
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% 6/26/13  dhb  Wrote it.

%% Default eye length
if (nargin < 5 || isempty(eyeLengthMm))
    eyeLengthMm = 17;
end
eyeLengthM = (10^-3)*eyeLengthMm;

%% Specify the limit (from table)
exposureDurationHours = stimulusDurationSecs/3600;
if (exposureDurationHours <= 2)
    limit_UWattsPerCm2 = 220;
else
    limit_UWattsPerCm2 = 220/(exposureDurationHours/2);
end

%% Convert radiance to retinal irradiance
%
% The standard says to do this with for a 7 mm pupil.  It does
% not give an eye length to assume.  We assume 17 mm.
pupilDiameterMm = 7;
pupilAreaMm2 = pi*((pupilDiameterMm/2)^2);
pupilAreaM2 = (10^-6)*pupilAreaMm2;

retIrradiance_WattsPerM2 = RadianceAndPupilAreaEyeLengthToRetIrradiance(radiance_WattsPerSrM2,S,pupilAreaM2,eyeLengthM);
retIrradiance_UWattsPerM2 = (10^6)*retIrradiance_WattsPerM2;
retIrradiance_UWattsPerCm2 = (10^-4)*retIrradiance_UWattsPerM2;

%% Get weighted sum.  The weighting function is zero outside the
% specified wavelength range, so we don't have to worry about the
% wavelength limits in the standard.  We do perform a sanity check
% that something got passed in the wavelength region of interest.
wls = SToWls(S);
index = find(wls >= 305 & wls <= 700, 1);
if (isempty(index))
    error('Should not call this routine with no spectral sampling between 305 and 700');
end
val_UWattsPerCm2 = sum(retIrradiance_UWattsPerCm2 .* weightingA);

