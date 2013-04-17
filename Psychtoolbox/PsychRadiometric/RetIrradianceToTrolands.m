function [trolands] =...
	RetIrradianceToTrolands(irradianceWatts, irradianceS, photopic, species, source)
% [trolands] =...
%  	RetIrradianceToTrolands(irradianceWatts,irradianceS,[photopic],[species],[source])
%
% Compute trolands from retinal irradiance in watts/um2-wlinterval.  The answer is
% returned in trolands/wlinterval and can be summed to get trolands.
%
% See Wyszecki and Stiles, 1982, p. 103 for the conversions.
% 
% Input variables: irradianceWatts - retinal irradiance in watts/um2-wlinterval.
%                  irradianceS - the wavelength sampling information for the relativeSpectrum.
%                  photopic - what kind of trolands: 'Photopic' (Default), 'JuddVos', 'Scotopic'. 
%                  species, source - passed directly to EyeLength to determine length of eye in mm.
%                     These values inherit the default behaviors of EyeLength.
%
% 07/18/03  dhb  Wrote it.
% 04/09/12  dhb  Debug.  This was apparently never quite finished.
%                Improve comments.

% Fill in default values
if (nargin < 3 || isempty(photopic))
	photopic = 'Photopic';
end
if (nargin < 4 || isempty(species))
	species = [];
end
if (nargin < 5 || isempty(source))
    source = [];
end
S = irradianceS;

% Convert spectrum to watts/deg2-wlinterval
eyeLengthMM = EyeLength(species,source);
mmPerDeg = DegreesToRetinalMM(1,eyeLengthMM);
mm2PerDeg2 = mmPerDeg^2;
irradianceMM2 = irradianceWatts*1e6;
irradianceDeg2 = irradianceMM2*mm2PerDeg2;

% Load appropriate V_lambda for phot/scot
switch (photopic)
	case 'Photopic'
		load T_xyz1931;
		T_vLambda = SplineCmf(S_xyz1931,T_xyz1931(2,:),S);
		clear T_xyz1931 S_xyz1931
		magicFactor = 2.242e12;
	case 'JuddVos'
		load T_xyzJuddVos;
		T_vLambda = SplineCmf(S_JuddVos,T_JuddVos(2,:),S);
		clear T_JuddVos S_JuddVos
		magicFactor = 2.242e12;
	case 'Scotopic'
		load T_rods;
		T_vLambda = SplineCmf(S_rods,T_rods,S);
		magicFactor = 5.581e12;
end

% Get trolands
trolands = T_vLambda*irradianceDeg2*magicFactor;
