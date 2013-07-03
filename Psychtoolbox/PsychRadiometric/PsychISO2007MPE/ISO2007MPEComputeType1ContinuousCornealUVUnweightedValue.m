function [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousCornealUVUnweightedValue(S,radiance_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2)
% [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousCornealUVUnweightedValue(S,radiance_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2)
%
% Compute the unweighted UV radiation for Type 1 instruments as given on page 7, Table 2, 
% 5.4.1.2.
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
% for microscopes and endoilluminators.
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% 6/25/13  dhb  Wrote it.

%% Specify the limit (from table)
%
% Limit in table is in mWatts, we convert to uWatts
% for consistency with the weighted routine.
exposureDurationHours = stimulusDurationSecs/3600;
if (exposureDurationHours <= 2)
    limit_UWattsPerCm2 = 1*(10^3);
else
    limit_UWattsPerCm2 = 1*(10^3)/(exposureDurationHours/2);
end

%% Convert radiance to corneal irradiance
cornealIrradiance_WattsPerM2 = RadianceAndDegrees2ToCornIrradiance(radiance_WattsPerSrM2,stimulusAreaDegrees2);
cornealIrradiance_UWattsPerM2 = (10^6)*cornealIrradiance_WattsPerM2;
cornealIrradiance_UWattsPerCm2 = (10^-4)*cornealIrradiance_UWattsPerM2;

%% Get sum betwen 360 and 400.
wls = SToWls(S);
index = find(wls >= 360 & wls <= 400);
if (isempty(index))
    error('Should not call this routine with no spectral sampling between 360 and 400');
end
val_UWattsPerCm2 = sum(cornealIrradiance_UWattsPerCm2(index));

