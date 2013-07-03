function cornealIrradiance_PowerPerArea = RadianceAndDistanceAreaToCornIrradiance(radiance_PowerPerSrArea,stimulusDistance,stimulusArea)
% cornealIrradiance_PowerPerArea = RadianceAndDistanceAreaToCornIrradiance(radiance_PowerPerSrArea,stimulusDistance,stimulusArea)
%
% Convert the radiance of a stimulus to corneal irradiance, given that we know the distance to the stimulus and the area
% of the stimulus.  Light power can be in your favorite units (Watts, quanta/sec) as can distance (m, cm, mm).  Area
% needs to be in units that are the square of your distance units, both for the radiance passed and the stimulus area
% passed. So, if radiance is in Watts/[cm2-sr] then distance needs to be in cm and irradiance will be in Watts/cm2.
%
% See also: RadianceAndDegrees2ToCornIrradiance, CornIrradianceAndDegrees2ToRadiance
%
% 2/20/13  dhb  Wrote it.

% Get total power coming off the stimulus.  According to Wyszecki and Stiles,
% Table 1(1.1), p. 2 this is called the radiant intensity.
radiantIntensity_PowerPerSr = radiance_PowerPerSrArea * stimulusArea;

% Figure out how much power per unit area by the time it arrives at the cornea.  To
% see that this is the right formula, we use the fact that 1 sr is the area
% subtended by the square of the radius of a sphere.  Here the radius is the
% distance between the emitting surface and the cornea, so that 1 sr is given
% by (stimulusDistance^2).  The radiant intensity is the power passing through 1
% sr at the cornea, and thus also the power passing through an area of (stimulusDistance^2).
% Dividing by this area gives the power per area.
cornealIrradiance_PowerPerArea = radiantIntensity_PowerPerSr/(stimulusDistance^2);

