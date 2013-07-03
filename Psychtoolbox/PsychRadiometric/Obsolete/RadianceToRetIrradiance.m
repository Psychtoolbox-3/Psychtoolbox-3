function irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,radianceS,pupilAreaMm2,eyeLengthMm)
% irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,radianceS,pupilAreaMm2,eyeLengthMm)
%
% Perform the geometric calculations necessary to convert a measurement of source
% radiance to corresponding retinal irradiance. 
%
%   Input radianceWattsPerM2Sr should be in units of power/m^2-sr-wlinterval.
%   Input radianceS gives the wavelength sampling information.
%   Input pupilAreaMm2 should be in units of mm^2.
%   Input eyeLengthMm should be the length of the eye in mm.
%   Output irradianceWattsPerUm2 is in units of power/um^2-wlinterval.
%
%   Light power may be expressed in watts or quanta-sec or in your
%   favorite units.  Indeed, it may also be passed as energy rather
%   than power.  
%
% This conversion does not take absorption in the eye into account,
% as this is more conveniently foldeded into the spectral absorbtance.
%
% See also: RadianceAndPupilAreaEyeLengthToRetIrradiance, radPupilAreaFromLum, EyeLength, RetIrradianceToRadiance.
%
% Note: This routine is now obsolete, as it mixes radiometric and unit conversions.  Preferred is to
% use RadianceAndPupilAreaEyeLengthToRetIrradiance and then take charge of your units in your
% calling code.
%
% 7/10/03  dhb  Wrote it.
% 11/06/03 dhb  Fixed comments about units, as per Lu Yin email.
% 3/29/12  dhb  Comment on output units was wrong.  It said power/um^2-sec-wlinterval
%               but the 'sec' part makes no sense given the 'power' in the numerator.
% 2/28/13  dhb  Make units clear in variable names.
% 3/6/13   dhb  Rewrite to use new conversion function.  Move to Obsolete directory.
%          dhb  Also improved variable naming.

%% Make sure we didn't break anything?
CHECKAGAINSTOLDCODE = 1;

%% Convert radiance units to match eye length and pupil units
radianceWattsPerMm2Sr = radianceWattsPerM2Sr*1e-6;

%% Convert to ret irradiance
irradianceWattsPerMm2 = RadianceAndPupilAreaEyeLengthToRetIrradiance(radianceWattsPerMm2Sr,radianceS,pupilAreaMm2,eyeLengthMm);

%% Convert irradiance units to uM2.
irradianceWattsPerUm2 = irradianceWattsPerMm2*1e-6;

%% Check what we get now against the original implementation
if (CHECKAGAINSTOLDCODE)
    
    % Convert spectral units to power/sr-mm^2-wlinterval
    radianceWattsPerMm2SrCheck = radianceWattsPerM2Sr*1e-6;
    
    % Define factor to convert radiance spectrum to retinal irradiance in watts/mm^2-wlinterval.
    % Commented out code shows the logic, which is short circuited by actual code.
    % but is conceptually convenient for doing the calculation.
    %  distanceToSourceMm = 100;
    %  fractionfSphere = pupilAreaMm2/4*pi*distanceToSourceMm^2;
    %  pupilAreaSR = 4*pi*fractionOfSphere;
    %  sourceAreaMm = (distanceToSourceMm^2)/(eyeLengthMm^2);
    %  conversionFactor = pupilAreaSR*sourceAreaMm;
    conversionFactor = pupilAreaMm2/(eyeLengthMm^2);
    irradianceWattsPerMm2Check = conversionFactor*radianceWattsPerMm2SrCheck;
    
    % Convert units to um^2 from mm^2 base.
    irradianceWattsPerUm2Check = irradianceWattsPerMm2Check*1e-6;
    
    % Check
    if (abs(irradianceWattsPerUm2 - irradianceWattsPerUm2Check) > 1e-16)
        error('New and old ways of computing this quantity do not agree.  Oops!');
    end
end

