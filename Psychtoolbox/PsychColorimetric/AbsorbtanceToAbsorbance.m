function [absorbanceSpectra, absorbanceSpectraWls] =...
    AbsorbtanceToAbsorbance(absorbtanceSpectra, absorbtanceSpectraWls, axialOpticalDensities, NORMALIZE)
% [absorbanceSpectra, absorbanceSpectraWls] =...
%   AbsorbtanceToAbsorbance(absorbtanceSpectra, absorbtanceSpectraWls, axialOpticalDensities, [NORMALIZE])
%
% This code inverts AbsorbanceToAbsorbtance.  You might want to do this if you were trying
% to take cone fundamentals and back down all the way to the component parts, so that you
% could for example vary the axial density and recompute the fundamentals.
%
% The absorbance/absorbtance terminology is described at the
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
% NORMALIZE (default true) causes this routine to normalize the returned absorbances to
% have a maximum of 1.
%
% Note, we now have ways of building up most fundamentals that we care about
% from constituant parts, and thus probably don't need to do that.  See
%   CIEConeConeFundamentlsTest, ComputeCIEConeFundamentals, DefaultPhotoreceptors, FillInPhotoreceptors,
%   IsomerizationsInEyeDemo.
%
% Originally written by HH, Copyright HH, Vista Lab, 2010
%
% 8/11/13  dhb  Moved into PTB, modified comments so as not to refer to non-PTB routines.
%          dhb  That this actually inverts is tested in IsomerizationsInEyeDemo.

% Some arg checks
if ~exist('absorbtanceSpectra','var'); help AbsorbtanceToAbsorbance; return; end
if ~exist('axialOpticalDensities','var'); disp('axialOpticalDensities is required.'); return; end
if ~exist('absorbtanceSpectraWls','var'); absorbtanceSpectraWls = []; end
if ~exist('NORMALIZE','var'); NORMALIZE = true; end

% Convert each entry
for i = 1:size(absorbtanceSpectra,1)
    % Normalize absorbtanceSpectance so that returned absorbance has peak of 1, no matter what
    % normalization was applied to the absorbtance spectrum.  
    if (NORMALIZE)
        absorbtanceSpectra(i,:) = absorbtanceSpectra(i,:) ./ max(absorbtanceSpectra(i,:)) * (1-10^-axialOpticalDensities(i));
    end
    
    % Invert the absorbance to absorbtance computation
    absorbanceSpectra(i,:) = (- 1 ./ axialOpticalDensities(i)) .* log10(absorbtanceSpectra(i,:) - 1);
end

% Deal with some possible numerical error
absorbanceSpectra = real(absorbanceSpectra); % remove small value of imaginary number

% Pass wavelengths back through
absorbanceSpectraWls = absorbtanceSpectraWls;

