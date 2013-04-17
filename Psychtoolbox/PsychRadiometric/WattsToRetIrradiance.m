function [irradiance, irradianceS] =...
	WattsToRetIrradiance(relativeSpectrum, relativeSpectrumS, readingInWatts, radiometer)
% [irradiance, irradianceS] =...
%		WattsToRetIrradiance(relativeSpectrum, relativeSpectrumS, readingInWatts, [radiometer])
%
% The assumption underlying this routine is that the relative spectrum of a light
% is available, as well as the total power of the light passing through an aperture
% of known size.  This is the situation in the apparatus we use in Sterling's lab.
%	The routine computes the irradiance (watts/um^2-wlinterval) from the relative spectrum
% (relative power, not relative quanta) and a measurement of the total power in watts,
% given the radiometer properties specified. 
%
% Radiometers do not typically have an ideal flat spectral response, and the actual
% spectral response is often provided as part of the instrument calibration.  In addition,
% the instruments are typically calibrated to provide an accurate reading for one
% particular reference spectrum, which is also specified.  This routine is set up
% to use this calibration information together with the known relative spectrum to
% provide a correct result.  The radiometer properties may be passed in the radiometer
% structure.  A default structure is set up that describes the radiometer used in
% the Sterling lab.  This structure also describes the collecting area of the measurement
% configuration.
%
% This routine could be easily modified to deal with a photometric head. 
% 
% Input variables: relativeSpectrum is the relative power as a function of wavelength.
%                  relativeSpectrumS is the wavelength sampling information for the relativeSpectrum.
%								   readingInWatts is the total power measured.
%                  radiometer is a structure describing the radiometer. (default = Sterling Lab's IL1400A).
%
% 04/29/03	lyin Wrote it with advice from DHB
% 05/06/03	lyin Put wls into variable: lightSource
% 05/06/03	lyin Correction for wavelength sampling interval
% 06/12/03	lyin Change the way, variable being passed
% 6/26/03   dhb  Change some names, also compute power per wavelength interval, not per nm.
% 6/30/03   dhb  Lots more changes.
% 7/08/03   dhb  Monochromatic ref spectrum default, as per manual.

% Fill in default radiometer properties.
if (exist('radiometer', 'var') ~= 1 || isempty(radiometer))
	% Structure containing radiometer calibration information.
	% The radiometer is calibrated to give the correct
	% reading in watts when a source with the relative spectrum
	% refspectrum is used, and that the intrinsic spectral efficiency
	% of the radiometer is the specified efficiency.
	%
	% The data below describe the Sterling Lab's IL1400A, whose reference
	% spectrum is a HeNe 633 nm line.	
	radiometer.pinholeDiameter = 1000;														    % um
	radiometer.efficiencyS = [350 10 46]; 														% nm
	radiometer.efficiency = [0.264 0.376 0.486 0.584 0.700...
													0.766 0.804 0.820 0.834 0.864...
													0.856 0.848 0.830 0.820 0.820...
													0.820 0.820 0.814 0.804 0.804...
													0.818 0.856 0.866 0.856 0.820...
													0.820 0.822 0.822 0.798 0.784...
													0.790 0.804 0.850 0.880 0.910...
													0.920 0.930 0.940 0.930 0.930...
													0.922 0.930 0.932 0.940 0.944 0.948];
	radiometer.refSpectrumS = [380 1 401];
	radiometer.refSpectrum = MakeMonoPrimary(633,radiometer.refSpectrumS);
	radiometer.measurementArea = pi * (radiometer.pinholeDiameter/2)^2;	% um^2
end

% Spline to common wavelength representation for computations
computeS = [380 1 401];
relativeSpectrum = SplineSpd(relativeSpectrumS,relativeSpectrum,computeS);
radiometer.efficiency = SplineCmf(radiometer.efficiencyS,radiometer.efficiency,computeS);
radiometer.refSpectrum = SplineSpd(radiometer.refSpectrumS,radiometer.refSpectrum,computeS);

% Figure out the internal calibration factor for the radiometer.  This is obtained by
% making sure that the refSpectrum, weighted by the efficiency, comes out at the
% correct reading.  Note that this factor only depends on the relative refSpectrum.
refSpectrumWatts = sum(radiometer.refSpectrum);
internalScaleFactor = refSpectrumWatts / (radiometer.efficiency*radiometer.refSpectrum);

% Solve for putting our spectrum into watts/wlinterval.
scaleFactor = readingInWatts / ( radiometer.efficiency*relativeSpectrum*internalScaleFactor );

% Get absolute spectrum in watts/wlinterval
irradiance2 = scaleFactor * relativeSpectrum;

% And now in watts/um^2-wlinterval
irradiance1 = irradiance2 / radiometer.measurementArea;

% Set returned wavelength sampling to match input.
irradiance = SplineSpd(computeS,irradiance1,relativeSpectrumS);
irradianceS = relativeSpectrumS;



