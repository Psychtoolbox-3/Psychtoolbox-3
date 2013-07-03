function [IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S_in,radiancein_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2,eyeLengthMm)
% [IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S_in,radiancein_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2,eyeLengthMm)
%
% Run all the checks that apply to a radiance measurement for the ISO 2007 MPE standard, for a Type 1 instrument and continuous exposure.
% Does not do the convergent beam check -- see comment in Contents.m about when we think this applies.
%
% Returns a flag set to 1 if the input exceeds any of the limits, and a structure with the computed ISO2007MPE
% values and corresponding limits.
%
% We use 17 mm eye length by default.  You can override the eye length by passing a different length in mm.
% This value is used in the conversion from radiance to retinal irradiance.
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% 6/28/13  dhb  Wrote it.

%% Eye length default
if (nargin < 5 || isempty(eyeLengthMm))
    eyeLengthMm = 17;
end

%% Wavelength sampling.
%
% Easiest to cover the whole range covered by the standard,
% that way we don't have to think.
S = [200 1 1301];

%% Read in the spectral functions used by the standard
[~,weightingR,weightingA,weightingS,wls_R,rawWeightingR,wls_A,rawWeightingA,wls_S,rawWeightingS] = ISO2007MPEGetWeighings(S);

%% Spline the input spectrum
radiance_WattsPerSrM2 = SplineSpd(S_in,radiancein_WattsPerSrM2,S);

%% Initialize limit flag
IsOverLimit = 0;

%% Corenal irradiance weighted UV limit
[ISO2007MPEStruct.cornealUVWeightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealUVWeightedLimit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousCornealUVWeightedValue(...
    S,radiance_WattsPerSrM2,weightingS,stimulusDurationSecs,stimulusAreaDegrees2);
if (ISO2007MPEStruct.cornealUVWeightedVal_UWattsPerCm2 >= ISO2007MPEStruct.cornealUVWeightedLimit_UWattsPerCm2)
    IsOverLimit = 1;
end

%% Corenal irradiance uweighted UV limit
[ISO2007MPEStruct.cornealUVUnweightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealUVUnweightedLimit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousCornealUVUnweightedValue(...
S,radiance_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2);
if (ISO2007MPEStruct.cornealUVUnweightedVal_UWattsPerCm2 >= ISO2007MPEStruct.cornealUVUnweightedLimit_UWattsPerCm2)
    IsOverLimit = 1;
end

%% Retinal irradiance weighted aphakic limit
[ISO2007MPEStruct.retIrradiancePCWeightedVal_UWattsPerCm2,ISO2007MPEStruct.retIrradiancePCWeightedLimit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousRetIrradiancePCWeightedValue(...
    S,radiance_WattsPerSrM2,weightingA,stimulusDurationSecs,eyeLengthMm);
if (ISO2007MPEStruct.retIrradiancePCWeightedVal_UWattsPerCm2 >= ISO2007MPEStruct.retIrradiancePCWeightedLimit_UWattsPerCm2)
    IsOverLimit = 1;
end

%% Radiance weighted aphakic limit
[ISO2007MPEStruct.radiancePCWeightedVal_UWattsPerSrCm2,ISO2007MPEStruct.radiancePCWeightedLimit_UWattsPerSrCm2] = ISO2007MPEComputeType1ContinuousRadiancePCWeightedValue(...
    S,radiance_WattsPerSrM2,weightingA,stimulusDurationSecs);
if (ISO2007MPEStruct.radiancePCWeightedVal_UWattsPerSrCm2 >= ISO2007MPEStruct.radiancePCWeightedLimit_UWattsPerSrCm2)
    IsOverLimit = 1;
end

%% Corneal irradiance unweighted IR limit
[ISO2007MPEStruct.cornealIRUnweightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealIRUnweightedLimit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousCornealIRUnweightedValue(...
    S,radiance_WattsPerSrM2,stimulusDurationSecs,stimulusAreaDegrees2);
if (ISO2007MPEStruct.cornealIRUnweightedVal_UWattsPerCm2 >= ISO2007MPEStruct.cornealIRUnweightedLimit_UWattsPerCm2)
    IsOverLimit = 1;
end

%% Retinal irradiance weighted thermal limit
[ISO2007MPEStruct.retIrradianceTHWeightedVal_UWattsPerCm2,ISO2007MPEStruct.retIrradianceTHWeightedLimit_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousRetIrradianceTHWeightedValue(...
    S,radiance_WattsPerSrM2,weightingR,stimulusDurationSecs,eyeLengthMm);
if (ISO2007MPEStruct.retIrradianceTHWeightedVal_UWattsPerCm2 >= ISO2007MPEStruct.retIrradianceTHWeightedLimit_UWattsPerCm2)
    IsOverLimit = 1;
end

%% Radiance weighted thermal limit
[ISO2007MPEStruct.radianceTHWeightedVal_UWattsPerSrCm2,ISO2007MPEStruct.radianceTHWeightedLimit_UWattsPerSrCm2] = ISO2007MPEComputeType1ContinuousRadianceTHWeightedValue(...
    S,radiance_WattsPerSrM2,weightingR,stimulusDurationSecs);
if (ISO2007MPEStruct.radianceTHWeightedVal_UWattsPerSrCm2 >= ISO2007MPEStruct.radianceTHWeightedLimit_UWattsPerSrCm2)
    IsOverLimit = 1;
end