function ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct)
% ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct)
%
% Useful printout of ISO 2007 analysis.  Formatted for pasting into Doku Wiki.
%
% 6/28/13  dhb  Wrote it.

if (IsOverLimit == 1)
    fprintf('  * Light is OVER ISO 2007 MPE Type 1 continuous limits\n');
else
    fprintf('  * Light is UNDER ISO 2007 MPE Type 1 continuous limits\n');
end
    
%% Corenal irradiance weighted UV limit
fprintf('  * Type 1 continuous corneal irradiance UV weighted (5.4.1.1)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/cm2)\n',ISO2007MPEStruct.cornealUVWeightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealUVWeightedLimit_UWattsPerCm2);

%% Corenal irradiance uweighted UV limit
fprintf('  * Type 1 continuous corneal irradiance UV unweighted (5.4.1.2)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/cm2)\n',ISO2007MPEStruct.cornealUVUnweightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealUVUnweightedLimit_UWattsPerCm2);

%% Retinal irradiance weighted aphakic limit
fprintf('  * Type 1 continuous aphakic retinal illumiance weighted (5.4.1.3.a)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/cm2)\n',ISO2007MPEStruct.retIrradiancePCWeightedVal_UWattsPerCm2,ISO2007MPEStruct.retIrradiancePCWeightedLimit_UWattsPerCm2);

%% Radiance weighted aphakic limit
fprintf('  * Type 1 continuous aphakic radiance weighted (5.4.1.3.b)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/[sr-cm2])\n',ISO2007MPEStruct.radiancePCWeightedVal_UWattsPerSrCm2,ISO2007MPEStruct.radiancePCWeightedLimit_UWattsPerSrCm2);

%% Corneal irradiance unweighted IR limit
fprintf('  * Type 1 continuous corneal irradiance IR unweighted (5.4.1.3.b)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/[sr-cm2])\n',ISO2007MPEStruct.cornealIRUnweightedVal_UWattsPerCm2,ISO2007MPEStruct.cornealIRUnweightedLimit_UWattsPerCm2);

%% Retinal irradiance weighted thermal limit
fprintf('  * Type 1 continuous thermal retinal illumiance weighted (5.4.1.3.a)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/cm2)\n',ISO2007MPEStruct.retIrradianceTHWeightedVal_UWattsPerCm2,ISO2007MPEStruct.retIrradianceTHWeightedLimit_UWattsPerCm2);

%% Radiance weighted thermal limit
fprintf('  * Type 1 continuous thermal radiance weighted (5.4.1.3.b)\n');
fprintf('    * Value: %0.3f, limit %0.3f (uWatts/[sr-cm2])\n',ISO2007MPEStruct.radianceTHWeightedVal_UWattsPerSrCm2,ISO2007MPEStruct.radianceTHWeightedLimit_UWattsPerSrCm2);