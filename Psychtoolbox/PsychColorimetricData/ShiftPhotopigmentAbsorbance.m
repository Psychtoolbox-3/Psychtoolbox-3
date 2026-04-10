function absorbance = ShiftPhotopigmentAbsorbance(S,absorbance,lambdaMaxShift,shiftMode)
% absorbance = ShiftPhotopigmentAbsorbance(S,absorbance,lambdaMaxShift,shiftMode)
%
% Function to shift photopigment absorbances.  Probably a reasonable
% approximation to biological reality for small shifts of lambda max.
%
% There are two ways of shifting photopigment absorbances:
% 'linear' - default
%   Shifts the photopigment absorbance on a linear wavelength axis.  This
%   is what Asano, Fairchild, & Bonde (2016), PLOS One, doi:
%   10.1371/journal.pone.0145671 do.
%
% 'log'
%   Shift the photopigment absorbance along a log wavelength axis to give it
%   a new lambda max.  This is the way Lamb (1995) suggested would lead to a
%   close approximation of biological shifts in photopigment lambda-max.
%
% Linearly extrapolates as necessary to preserve wavelength sampling of
% input.
%
% 2/8/16  dhb, ms  Wrote it.
% 2/10/16 dhb, ms  Do shifting on a fine (0.25 nm) wavelength spacing.
% 2/12/16 ms       Implement linear and log shifting.

% Assume linear shifting as default
if (nargin < 4 | isempty(shiftMode))
    shiftMode = 'linear';
end

% Check
if (length(lambdaMaxShift) ~= size(absorbance,1))
    error('Number of absorbances and length of shift vector do not match');
end

% Do the shifting on a fine wavelength spacing scale, so as to avoid
% unexpected wavelength quantization artifacts as a result of the shift.
wls = MakeItWls(S);
fineSpacing = 0.25;
if (S(2) < fineSpacing)
    error('Not a good idea to run this function with an input wavelength spacing less than 0.25 nm');
end
wlsFine = (wls(1):fineSpacing:wls(end))';
SFine = WlsToS(wlsFine);
absorbanceFine = SplineCmf(S,absorbance,SFine);

switch shiftMode
    case 'linear'
        for ii = 1:size(absorbanceFine,1);
            absorbanceFine(ii, :) = interp1(wlsFine+lambdaMaxShift(ii), absorbanceFine(ii, :), wlsFine, 'linear', 'extrap');
        end
    case 'log'
        % Normalize the wave number and log wavelength
        [~, maxIdx] = max(absorbanceFine, [], 2);
        for ii = 1:size(absorbanceFine,1);
            logWavelengthsNorm = log10(wlsFine) - log10(wlsFine(maxIdx(ii)));
            logWavelengthsNew = logWavelengthsNorm + log10(wlsFine(maxIdx(ii))+lambdaMaxShift(ii));
            wlsFineNew = 10.^logWavelengthsNew;
            absorbanceFine(ii, :) = interp1(wlsFineNew, absorbanceFine(ii, :), wlsFine, 'linear','extrap');
        end
end

% Spline the shifted function back to the original wavelength spacing
absorbance = SplineCmf(SFine,absorbanceFine,S);