function [absorbtanceSpectra, absorbtanceSpectraWls] =...
	AbsorbanceToAbsorbtance(absorbanceSpectra, absorbanceSpectraWls, axialOpticalDensities)
% [absorbtanceSpectra, absorbtanceSpectraWls] =...
%   AbsorbanceToAbsorbtance(absorbanceSpectra, absorbanceSpectraWls, axialOpticalDensities)
%
% Convert pigment absorbance spectra into absorbtance spectra, using the peak axial
% optical density.  The absorbance/absorbtance terminology is described at the
% CVRL web page, http://cvrl.ucl.ac.uk.  Wyszecki and Stiles refere to absorbance
% the absorption coefficient (p. 588).
%
% Both absorbtance spectra and absorbance spectra describe quantal absorption.
%
% Absorbance spectra are normalized to a peak value of 1.
% Absorbtance spectra are the proportion of quanta actually absorbed.
%
% Equation: absorbtanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
%
% Multiple spectra may be passed in the rows of absorbanceSpectra.  If
% so, then the same number of densities should be passed in the vector
% axialOpticalDensities, and multiple answers are returned in the rows
% of absorbtanceSpectra.
%
% Wavelength information may be in any of the available Psychtoolbox representations,
% and the returned wavelength information is in the same format as passed.
%
% 04/29/03 lyin 	Wrote wrote with advice from dhb
% 04/30/03 lyin 	Reorganize the variable
% 06/12/03 lyin 	Change the way variable being passed
% 06/23/03 dhb		Check dimensions of spectra and density.
% 06/30/03 dhb      Change to toolbox convention, put sensitivity like stuff in rows.
% 08/11/13 dhb      Fix comment to reflect row convention change made in 2003.  Slowly but surely we fix things up.

% Check that dimensions match properly
if (size(absorbanceSpectra,1) ~= length(axialOpticalDensities))
	error('Number of spectra does not match number of densities');
end

% Equation: absorbtanceSpectra = 1 - 10.^(-OD * absorbanceSpectra)
absorbtanceSpectra = 1 - 10.^(-diag(axialOpticalDensities)*absorbanceSpectra);

% Wls of absorbtanceSpectra
absorbtanceSpectraWls = absorbanceSpectraWls;
