function retIrradiance_PowerPerArea = RadianceAndPupilAreaEyeLengthToRetIrradiance(radiance_PowerPerAreaSr,radianceS,pupilArea,eyeLength)
% retIrradiance_PowerPerArea = RadianceAndPupilAreaEyeLengthToRetIrradiance(radiance_PowerPerAreaSr,radianceS,pupilArea,eyeLength)
%
% Perform the geometric calculations necessary to convert a measurement of source
% radiance to corresponding retinal irradiance. 
%
% Let x be the units of distance (m, cm, mm, um, etc.)
%
%   Input radiance_PowerPerAreaSr should be in units of power/x^2-sr-wlinterval.
%   Input radianceS gives the wavelength sampling information (nm).
%   Input pupilArea should be in units of x^2.
%   Input eyeLength should be the length of the eye in units of x
%
%   Output retIrradiance_PowerPerArea is in units of power/x^2-wlinterval.
%
%   Light power may be expressed in watts or quanta-sec or in your
%   favorite units.  Indeed, it may also be passed as energy rather
%   than power.  
%
% This conversion does not take absorption in the eye into account,
% as this is more conveniently foldeded into the spectral absorbtance.
%
% See also: PsychRadiometric, RetIrradianceAndPupilAreaEyeLengthToRadiance, PupilAreaFromLum, EyeLength.
%
% 3/6/13  dhb  Wrote it.

% Define factor to convert radiance spectrum to retinal irradiance
% Commented out code shows the logic, which is short circuited by actual code.
% but is conceptually convenient for doing the calculation.
%  distanceToSource = 100;
%  fractionfSphere = pupilArea/4*pi*distanceToSource^2;
%  pupilAreaSR = 4*pi*fractionOfSphere;
%  sourceArea = (distanceToSource^2)/(eyeLength^2);
%  conversionFactor = pupilAreaSR*sourceArea;
conversionFactor = pupilArea/(eyeLength^2);
retIrradiance_PowerPerArea = conversionFactor*radiance_PowerPerAreaSr;
