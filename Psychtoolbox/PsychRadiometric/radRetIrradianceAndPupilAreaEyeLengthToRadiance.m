function radiance_PowerPerAreaSr = radRetIrradianceAndPupilAreaEyeLengthToRadiance(irradiance_PowerPerArea,irradianceS,pupilArea,eyeLength)
% radiance_PowerPerAreaSr = radRetIrradianceAndPupilAreaEyeLengthToRadiance(irradiance_PowerPerArea,irradianceS,pupilArea,eyeLength)
%
% Perform the geometric calculations necessary to convert a measurement of retinal
% irradiance to the source radiance that would produce it.
%
% Perform the geometric calculations necessary to convert a measurement of source
% radiance to corresponding retinal irradiance. 
%
% Let x be the units of distance (m, cm, mm, um, etc.)
%
%   Input irradiance_PowerPerArea is in units of power/x^2-wlinterval.
%   Input irradianceS gives the wavelength sampling information.
%   Input pupilArea should be in units of x^2.
%   Input eyeLength should be the length of the eye in x.
%
%   Output radiance_PowerPerAreaSr is in units of power/x^2-sr-wlinterval.
%
%   Light power may be expressed in watts or quanta-sec or in your
%   favorite units.  Indeed, it may also be passed as energy rather
%   than power.  
%
% This conversion does not take absorption in the eye into account,
% as this is more conveniently foldeded into the spectral absorbtance.
%
% See also: PsychRadiometric, radRadianceAndPupilAreaEyeLengthToRetIrradiance, PupilAreaFromLum, EyeLength.
%
% 3/6/13  dhb  Wrote it.


% Define factor to convert radiance spectrum to retinal irradiance
% and apply this in the opposite direction.  See
% radRadianceAndPupilAreaEyeLengthToRetIrradiance for the
% conversion logic.
conversionFactor = pupilArea/(eyeLength^2);
radiance_PowerPerAreaSr = irradiance_PowerPerArea/conversionFactor;
