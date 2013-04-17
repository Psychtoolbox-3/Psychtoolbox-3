function radianceWattsPerM2Sr = RetIrradianceToRadiance(irradianceWattsPerUm2,irradianceS,pupilAreaMm2,eyeLengthMm)
% radianceWattsPerM2Sr = RetIrradianceToRadiance(irradianceWattsPerUm2,irradianceS,pupilAreaMm2,eyeLengthMm)
%
% Perform the geometric calculations necessary to convert a measurement of retinal
% irradiance to the source radiance that would produce it.
%
%   Input irradianceWattsPerUm2 is in units of power/um^2-wlinterval.
%   Input irradianceS gives the wavelength sampling information.
%   Input pupilAreaMm2 should be in units of mm^2.
%   Input eyeLengthMm should be the length of the eye in mm.
%   Output radianceWattsPerM2Sr is in units of power/m^2-sr-wlinterval.
%
%   Light power may be expressed in watts or quanta-sec or in your
%   favorite units.  Indeed, it may also be passed as energy rather
%   than power.  
%
% This conversion does not take absorption in the eye into account,
% as this is more conveniently foldeded into the spectral absorbtance.
%
% See also: radRetIrradianceAndPupilAreaEyeLengthToRadiance
%
% Note: This routine is now obsolete, as it mixes radiometric and unit conversions.  Preferred is to
% use radRadianceAndPupilAreaEyeLengthToRetIrradiance and then take charge of your units in your
% calling code.
%
% 2/28/13  dhb  Wrote it.
% 3/6/13   dhb  Rewrite to use new conversion function.  Move to Obsolete directory.
%          dhb  Also improved variable naming.

%% Make sure we didn't break anything?
CHECKAGAINSTOLDCODE = 1;

%% Convert units from um^2 to mm^2 base.
irradianceWattsPerMm2 = irradianceWattsPerUm2*1e6;

%% Convert to radiance
radianceWattsPerMm2Sr = radRetIrradianceAndPupilAreaEyeLengthToRadiance(irradianceWattsPerMm2,irradianceS,pupilAreaMm2,eyeLengthMm);

%% Convert radiance to power/sr-m^2-wlinterval
radianceWattsPerM2Sr = radianceWattsPerMm2Sr*1e6;

%% Check what we get now against the original implementation
if (CHECKAGAINSTOLDCODE)
    % Convert units from um^2 to mm^2 base.
    irradianceWattsPerMm2Check = irradianceWattsPerUm2*1e6;
    
    % Define factor to convert radiance spectrum to retinal irradiance in watts/mm^2-wlinterval.
    % We then apply this in the opposite direction.  See RadianceToRetIrradiance for the
    % conversion logic.
    conversionFactor = pupilAreaMm2/(eyeLengthMm^2);
    radianceWattsPerMm2SrCheck = irradianceWattsPerMm2Check/conversionFactor;
    
    % Convert radiance to power/sr-m^2-wlinterval
    radianceWattsPerM2SrCheck = radianceWattsPerMm2SrCheck*1e6;
    
    % Check
    if (abs(radianceWattsPerM2Sr - radianceWattsPerM2SrCheck) > 1e-16)
        error('New and old ways of computing this quantity do not agree.  Oops!');
    end
end


