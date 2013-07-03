% ISO2007MPEBasicTest
%
% ****************************************************************************
% IMPORTANT: Before using the ISO2007MPE routines, please see the notes on usage
% and responsibility in PsychISO2007MPE/Contents.m (type "help PsychISO2007MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Test code for our implementation of the ISO 2007 broadband MPE standard.
%
% We don't have any known test cases to check, particularly since we can
% only measure light in the visible.  Here we verify that a bright sunlight
% measured in Philly is below the MPE, but not by all that much.  This is
% broadly consistent with what we find when we compare the same light to the ANSI
% standard for laser light (with a few assumptions to apply that standard to
% broadband light.)
%
% 6/26/13  dhb  Wrote it.

%% Clear and close
clear; close all;

%% Wavelength sampling.
%
% Easiest to cover the whole range covered by the standard,
% that way we don't have to think.
S = [200 1 1301];

%% Load CIE functions.   
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);

%% Load in a test spectrum
%
% This is a bright sunlight measured through a window and off
% of a white piece of paper in Philly.
% We only have measurements between 380 and 780 nm.  In this
% example, we zero extend when we spline to the whole range.
load spd_phillybright
spd_phillybright = SplineSpd(S_phillybright,spd_phillybright,S,0);
photopicLuminancePhillyBrightCdM2 = T_xyz(2,:)*spd_phillybright;
PLOT_SPECTRUM = 0;
if (PLOT_SPECTRUM)
    figure;
    plot(SToWls(S),spd_phillybright,'r','LineWidth',2);
    xlabel('Wavelength (nm)');
    ylabel('Radiance (Watts/[sr-m2-wlinterval]');
end

%% Specify stimulus parameters
stimulusDiameterDegrees = 10;
stimulusAreaDegrees2 = pi*((stimulusDiameterDegrees/2)^2);
stimulusDurationSecs = 60*60;

%% Eye length, for conversion from radiance to retinal irradiance
%
% This will default to 17 if you don't pass it to the compute routine.
eyeLengthMm = 17;

%% Plot of weighting functions
[wls,weightingR,weightingA,weightingS,wls_R,rawWeightingR,wls_A,rawWeightingA,wls_S,rawWeightingS] = ISO2007MPEGetWeighings(S);
figure; clf; set(gcf,'Position',[400 500 1400 550]);
subplot(1,3,1); hold on
plot(wls_R,rawWeightingR,'r','LineWidth',3);
plot(wls,weightingR,'b:','LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('R_lambda')
title('Thermal Hazard Function');
xlim([200 1500]);
subplot(1,3,2); hold on
plot(wls_A,rawWeightingA,'r','LineWidth',3);
plot(wls,weightingA,'b:','LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('A_lambda')
title('Aphakic Photochemical Hazard Function');
xlim([200 1500]);
subplot(1,3,3); hold on
plot(wls_S,rawWeightingS,'r','LineWidth',3);
plot(wls,weightingS,'b:','LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('S_lambda')
title('UV Radiation Hazard Function');
xlim([200 1500]);

%% Do the computations
fprintf('\nRunning tests for a sunlight measured in Philadelphia\n\n');
[IsOverLimit,ISO2007MPEStruct] = ISO2007MPECheckType1ContinuousRadiance(S,spd_phillybright,stimulusDurationSecs,stimulusAreaDegrees2,eyeLengthMm);
ISO2007MPEPrintAnalysis(IsOverLimit,ISO2007MPEStruct)

%% Anterior segment limit for convergent beams
% This just makes sure the routine properly throws an error.
try
    fprintf('\n\nTesting error catch for anterior segment limit for convergent beams\n');
    [ISO2007MPEStruct.val8_UWattsPerCm2,limit8_UWattsPerCm2] = ISO2007MPEComputeType1ContinuousAntConvrgUnweightedValue(S,NaN,stimulusDurationSecs);   
catch
end