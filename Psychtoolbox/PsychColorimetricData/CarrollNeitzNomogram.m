function T_absorbance = CarrollNeitzNomogram(S, lambdaMax)
% T_absorbance = CarrollNeitzNomogram(S, lambdaMax)
%
% Compute normalized absorbance according to the nomogram of
% Carroll, McMahon, Neitz, & Neitz (2000), J. Opt. Soc. Am. A, 17(3), 499-509.
%
% T_absorbance contains the absorbance (not log absorbance).
% The CVRL convention is used: absorbance = log(I_incident/I_transmitted),
% so values are positive, with a peak of 1.  The peak is normalized
% explicitly to 1 when wavelength spacing S(2) <= 1 nm.
%
% The result is in quantal units.  To get sensitivity in energy units,
% apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of lambda max values, in which
% case T_absorbance is nLambdaMax x nWls.
%
% The nomogram is parameterized relative to a reference lambda max of
% 558.5 nm.  The internal shifted frequency variable is
%   x = lambdaMax / (558.5 * wl)
% which is analogous to the ratio form used by the Govardovskii nomogram.
%
% Reasonable lambda max values for (L, M, S): 557.5, 530, 420 nm.
% These values are from Neitz & Neitz (2011), who report good agreement
% between this template and CIE cone fundamentals derived from color matching
% at these peaks; see their Figure 3 and surrounding text.
%   Neitz, J. & Neitz, M. (2011). The genetics of normal and defective
%   color vision. Vision Research, 51, 633-651.
%
% See also CarrollNeitzNomogramTest, PhotpigmentNomogramDemo.

% History:
%   2026-04-10  dhb  Wrote it, based on spectsens.m from Neitz lab website.

% Wavelength sampling
S = MakeItS(S);
wls = MakeItWls(S);
nWls  = length(wls);
nT    = length(lambdaMax);
T_absorbance = zeros(nT, nWls);

% -----------------------------------------------------------------------
% Empirical parameters for the Carroll-Neitz nomogram.
% Names follow the role of each parameter rather than the original
% single-letter labels A-Y used in spectsens.m.
% -----------------------------------------------------------------------

% Primary tanh component (exTemp2 in spectsens.m)
tanhAmp = 0.417050601;       % overall amplitude
tanhCtr = 0.002072146;       % center in shifted-frequency units
tanhWid = 0.000163888;       % half-width

% Log-tanh component (exTemp1 in spectsens.m)
logTanhAmp = -16.05774461;   % amplitude (negative, so -E+E*tanh... is a rising sigmoid)
logTanhCtr =  0.001575426;   % center
logTanhWid =  5.11376e-05;   % half-width
logTanhOff = -1.922880605;   % additive offset in log space

% Gaussian correction terms (negative corrections, exTemp3-5)
gaussCtr1 = 0.00157981;   gaussSig1 = 6.58428e-05;  gaussAmp1 = 6.68402e-05;
gaussCtr2 = 0.002310442;  gaussSig2 = 7.31313e-05;  gaussAmp2 = 1.86269e-05;
gaussCtr3 = 0.002008124;  gaussSig3 = 5.40717e-05;  gaussAmp3 = 5.14736e-06;

% Gaussian correction terms (positive corrections, exTemp6-8)
gaussCtr4 = 0.001455413;  gaussSig4 = 4.21764e-05;  gaussAmp4 = 4.800e-06;
gaussCtr5 = 0.001809022;  gaussSig5 = 3.86677e-05;  gaussAmp5 = 2.990e-05;
gaussCtr6 = 0.001757315;  gaussSig6 = 1.47344e-05;  gaussAmp6 = 1.510e-05;

gaussNorm = 1/sqrt(2*pi);

% -----------------------------------------------------------------------
% Compute absorbance for each lambda max
% -----------------------------------------------------------------------
for ii = 1:nT

    % Shifted frequency variable (nWls x 1 column vector).
    % log10(1/wl) - [log10(1/lambdaMax) - log10(1/558.5)] simplifies to:
    x = lambdaMax(ii) ./ (558.5 .* wls);

    % Log absorbance: sum of two tanh terms and six Gaussian corrections
    logTanhTerm  = log10(-logTanhAmp + logTanhAmp*tanh(-((x) - logTanhCtr)/logTanhWid)) + logTanhOff;
    tanhTerm     = tanhAmp * tanh(-((x) - tanhCtr)/tanhWid);

    gaussTerm1   = -(gaussAmp1/gaussSig1) * gaussNorm .* exp(-0.5*((x - gaussCtr1)/gaussSig1).^2);
    gaussTerm2   = -(gaussAmp2/gaussSig2) * gaussNorm .* exp(-0.5*((x - gaussCtr2)/gaussSig2).^2);
    gaussTerm3   = -(gaussAmp3/gaussSig3) * gaussNorm .* exp(-0.5*((x - gaussCtr3)/gaussSig3).^2);
    gaussTerm4   =  (gaussAmp4/gaussSig4) * gaussNorm .* exp(-0.5*((x - gaussCtr4)/gaussSig4).^2);
    gaussTerm5   =  (gaussAmp5/gaussSig5) * gaussNorm .* exp(-0.5*((x - gaussCtr5)/gaussSig5).^2) / 10;
    gaussTerm6   =  (gaussAmp6/gaussSig6) * gaussNorm .* exp(-0.5*((x - gaussCtr6)/gaussSig6).^2) / 100;

    logAbsorbance = logTanhTerm + tanhTerm + ...
                    gaussTerm1 + gaussTerm2 + gaussTerm3 + ...
                    gaussTerm4 + gaussTerm5 + gaussTerm6;

    % Convert log absorbance to absorbance and store as a row
    T_absorbance(ii,:) = (10.^logAbsorbance)';

    % Normalize to peak 1 when wavelength spacing is fine enough that the
    % sampled peak reliably captures the true maximum.
    if (S(2) <= 1)
        T_absorbance(ii,:) = T_absorbance(ii,:) / max(T_absorbance(ii,:));
    end
end

