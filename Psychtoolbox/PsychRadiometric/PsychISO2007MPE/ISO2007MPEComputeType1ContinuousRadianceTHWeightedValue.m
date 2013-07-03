function [val_UWattsPerSrCm2,limit_UWattsPerSrCm2] = ISO2007MPEComputeType1ContinuousRadianceTHWeightedValue(...
    S,radiance_WattsPerSrM2,weightingR,stimulusDurationSecs)
%[val_UWattsPerSrCm2,limit_UWattsPerSrCm2] = ISO2007MPEComputeType1ContinuousRadianceTHWeightedValue(...
%    S,radiance_WattsPerSrM2,weightingR,stimulusDurationSecs)
%
% Compute the weighted thermal radiance for Type 1 instruments as given on page 9, Table 2, 
% 5.4.1.6.b.  Note.  The limit specified there is 6 W/[sr-cm2].  I reduced to 5.88 because
% this is the worst-case Type 2 limit (largest retinal diameter spot, see Table 4, 5.5.1.5b),
% and it seemed more conservative to do so.
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
% The indicates that radiance should be measured through a 7 mm aperture at the cornea.
% I believe this refers to the case where you are measuring radiance by measuring radiant
% flux through two apertures at a known distance apart.  It would then make sense that you'd
% want to know the radiance defined by the direction subtended by a 7 mm aperture at the
% cornea, e.g. right where a large pupil would be.  If you measure using some other device
% (e.g. a PhotoResearch PR-XXX that directly obtains radiance and do so from the eye position,
% that also seems reasonable.
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% 6/26/13  dhb  Wrote it.

%% Specify the limit (from table, see note above about why
% 5.88 rather than 6.)
exposureDurationHours = stimulusDurationSecs/3600;
if (exposureDurationHours <= 2)
    limit_UWattsPerSrCm2 = 5.88*(10^6);
else
    limit_UWattsPerSrCm2 = 5.88*(10^6)/(exposureDurationHours/2);
end

%% Unit conversion
radiance_UWattsPerSrM2 = (10^6)*radiance_WattsPerSrM2;
radiance_UWattsPerSrCm2 = (10^-4)*radiance_UWattsPerSrM2;

%% Get weighted sum.  The weighting function is zero outside the
% specified wavelength range, so we don't have to worry about the
% wavelength limits in the standard.  We do perform a sanity check
% that something got passed in the wavelength region of interest.
wls = SToWls(S);
index = find(wls >= 380 & wls <= 1400, 1);
if (isempty(index))
    error('Should not call this routine with no spectral sampling between 380 and 1400');
end
val_UWattsPerSrCm2 = sum(radiance_UWattsPerSrCm2 .* weightingR);

