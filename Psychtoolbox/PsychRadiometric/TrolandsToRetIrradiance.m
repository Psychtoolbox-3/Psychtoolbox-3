function [irradianceWattsPerUm2, irradianceS] =...
	TrolandsToRetIrradiance(relativeSpectrum, relativeSpectrumS, trolands, photopic, species, source)
% [radianceWattsPerM2Sr, irradianceS] =...
%  	TrolandsToRetIrradiance(relativeSpectrum, relativeSpectrumS, trolands, [photopic], [species], [source])
%
% The assumption underlying this routine is that the relative spectrum of a light
% is available, as well as the retinal illuminance in trolands.
%
% The routine computes the irradiance (watts/um^2-wlinterval) from the relative spectrum
% (relative power, not relative quanta) and the number of trolands.
%
% See Wyszecki and Stiles, 1982, p. 103 for the conversions.
% 
% Input variables: relativeSpectrum - the relative power as a function of wavelength.
%                  relativeSpectrumS - the wavelength sampling information for the relativeSpectrum.
%                  trolands - the number of trolands.
%                  photopic - what kind of trolands: 'Photopic' (Default), 'JuddVos', 'Scotopic'. 
%                  species - what species determins eye length: 'Human' (Default), 'Monkey'.
%                  source - source for eye length estimate, passed directly to EyeLength and inherits its default.
%
% 07/18/03  dhb         Wrote it.
% 1/26/04   ly, dhb     Fix JuddVos path through switch.
% 7/16/13   dhb         Comment and code cleaning, minor.

% Fill in default values
if (nargin < 4 || isempty(photopic))
	photopic = 'Photopic';
end
if (nargin < 5 || isempty(species))
	species = 'Human';
end
S = relativeSpectrumS;

% Load appropriate V_lambda for phot/scot
switch (photopic)
	case 'Photopic'
		load T_xyz1931;
		T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
		clear T_xyz1931 S_xyz1931
		magicFactor = 2.242e12;
	case 'JuddVos'
		load T_xyzJuddVos;
		T_vLambda = SplineCmf(S_xyzJuddVos,T_xyzJuddVos(2,:),S);
		clear T_JuddVos S_JuddVos
		magicFactor = 2.242e12;
	case 'Scotopic'
		load T_rods;
		T_vLambda = SplineCmf(S_rods,T_rods,S);
		magicFactor = 5.581e12;
end

% Convert relative spectrum into watts/deg^2 on retina.
% We know that trolands = k*(T_vLambda*relativeSpectrum)*magicFactor,
% were k is a factor that puts the relative spectrum into units of
% retinal irradiance in watts/deg^2-wlinterval
k = trolands/(T_vLambda*relativeSpectrum*magicFactor);
irradianceDeg2 = k*relativeSpectrum;

% Convert to units of um2
eyeLengthMM = EyeLength(species,source);
mmPerDeg = DegreesToRetinalMM(1,eyeLengthMM);
mm2PerDeg2 = mmPerDeg^2;
irradianceMM2 = irradianceDeg2/mm2PerDeg2;
irradianceWattsPerUm2 = irradianceMM2*1e-6;
irradianceS = S;
