% RadiometricConversionsTest
%
% Test out the radiometric conversion routines.
%
% This is not a complete test.  And, I thought I wrote such a thing
% once before but can't currently find it.
%
% 6/23/13  dhb  Wrote it.

%% Clear
clear; close all;

%% Define a monochromatic radiance
S = [500 1 1];
radianceWattsPerSrM2 = ones(S(3),1);

%% Define eye parameters
eyeLengthM = 0.017;
pupilDiameterM = 0.003;
pupilAreaM2 = 2*pi*((pupilDiameterM/2)^2);

%% Convert radiance to retinal irradiance
retIrradianceWattsPerM2 = RadianceAndPupilAreaEyeLengthToRetIrradiance(radianceWattsPerSrM2,S,pupilAreaM2,eyeLengthM);

%% Convert retinal irradiance to units of degrees2.
retIrradianceWattsPerDegrees2 = RetIrradiancePerAreaAndEyeLengthToRetIrradiancePerDegrees2(retIrradianceWattsPerM2,eyeLengthM);

%% Convert back to units of area and check that we get what we put in.
retIrradianceWattsPerM2Check = RetIrradiancePerDegrees2AndEyeLengthToRetIrradiancePerArea(retIrradianceWattsPerDegrees2,eyeLengthM);
fprintf('Retinal irradiance: %0.3g Watts/M2 (%0.3g check), %0.3g Watts/deg2\n', ...
    retIrradianceWattsPerM2,retIrradianceWattsPerM2Check,retIrradianceWattsPerDegrees2);
if (abs(retIrradianceWattsPerM2-retIrradianceWattsPerM2Check) > 1e-10)
    error('Retinal irradiance unit calculations do not invert');
else
    
end


