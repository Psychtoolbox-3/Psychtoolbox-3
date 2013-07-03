function [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousAntConvrgUnweightedValue(S,irradiance_UWattsPerCm2,stimulusDurationSecs)
% [val_UWattsPerCm2,limit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousAntConvrgUnweightedValue(S,irradiance_uWattsPerCm2,stimulusDurationSecs)
%
% Compute the unweighted radiation for anterior segment convergent beam value, for Type 1 instruments as given on page 8, Table 2, 
% 5.4.1.5.  This limit applies only to convergent beams, whatever they are.  I'm guessing this is what a Maxwellian view produces,
% however, as in that case the beam is brought to a waist inside the eye.  In this case, the limit applies to the irradiance
% at the beam waist over the 1 mm diameter aperture with the highest irradiance.
%
% Unlike all the other routines in this suite, the key quantity is probably not best passed as stimulus radiance.  Rather, for
% this type of optical system, the relevant irradiance should be measured directly.
%
% NOTE: This routine has not been tested, since I (DHB) don't have any instrument currently that produces the relevant type of
% stimulus.  I added it as a placeholder for convenience.  There is an error statement to keep it from being used
% until it is tested.
%
% Input spectrum is radiance in units of UWatts/[m2-wlinterval].
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
% 6/28/13  dhb  Wrote it.

%% Error statement to indicate draft status of this routine.
fprintf('\tISO2007MPEComputeType1ContinuousAntConvrgUnweightedValue\n');
fprintf('\t\tThis routine has been drafted but not tested.  It might\n');
fprintf('\t\tbe correct.  Or it might not.  Develop a reasonable\n');
fprintf('\t\ttest before using.\n');
error('Not yet ready for use');

%% Specify the limit (from table)
%
% Limit in table is in Watts, we convert to uWatts
% for consistency with the weighted routine.
exposureDurationHours = stimulusDurationSecs/3600;
if (exposureDurationHours <= 2)
    limit_UWattsPerCm2 = 4*(10^6);
else
    limit_UWattsPerCm2 = 4*(10^6)/(exposureDurationHours/2);
end

%% Get sum betwen 380 and 1200.
wls = SToWls(S);
index = find(wls >= 380 & wls <= 1200);
if (isempty(index))
    error('Should not call this routine with no spectral sampling between 380 and 1200');
end
val_UWattsPerCm2 = sum(irradiance_UWattsPerCm2(index));

