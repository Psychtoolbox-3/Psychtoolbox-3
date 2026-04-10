
function T_absorbance = StockmanSharpeNomogram(S,lambdaMax)
% T_absorbance = StockmanSharpeNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided by Stockman and Sharpe:
%   See Stockman & Sharpe (2000), p. 1730 or http://www.cvrl.org/.
%
% The lmax values of the fitted templates that best fits 
% to the Stockman and Sharpe (2000) S-, M- and L-cone photopigment
% spectra are 420.7, 530.3 and 558.9 nm for the
% S-, M- and L-cones, respectively;
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta()
% (not QuantaToEnergy, because here you are converting sensitivity
% rather than spectra.)
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% This routine converts the log10 absorbance computed by the nomogram
% formulae into absorbance.
%
% Note from DHB. By eye, this nomogram gives a good fit to the log10 LMS pigment
% absorbance spectra when you plot them on a log scale over 8 log units,
% as in Stockman & Sharpe (2000), Figure 12.  The fit does not look quite so
% good in my hands on a linear scale or an expanded log scale.
% I think the deviations occur in part because the tabulated
% photopigment absorbances for the L are a mixture of the ser/ala 
% pigments, and the nomogram was developed in part on the basis of fitting these
% as if they corresponded to a single lambda max value. It may also be that
% the template was built by minimizing the error in log sensitivity, and this
% would more heavily weight the long wavelength limbs of the pigments, where
% the linear sensitivity is essentially zero.  In any case, though, if you
% are working in the land of the CIE 170-1:2006 fundamentals, this is the
% probably the best current nomogram to use.
%
% See ComputeCIEConeFundamentals, CIEConeFundamentalsTest,
% FitConeFundamentalsFromNomogram, FitConeFundamentalsTest
%
% 5/8/99	dhb  Started writing it.
% 10/27/99	dhb  Added error return to prevent premature use of this routine.
% 7/18/03   dhb  Finished it off.
% 8/13/11   dhb  Improved comments.  Double check polynomial coefficients.

% Parameters
a = -188862.970810906644;
b = 90228.966712600282;
c = -2483.531554344362;
d = -6675.007923501414;
e = 1813.525992411163;
f = -215.177888526334;
g = 12.487558618387;
h = -0.289541500599;

% Set up and apply formula
S = MakeItS(S);
wls = MakeItWls(S)';
nWls = length(wls);
nT = length(lambdaMax);
T_absorbance = zeros(nT,nWls);
for i = 1:nT
	% Get lambda max
	theMax = lambdaMax(i);
	
	% Need to normalize wavelengths
	logWlsNorm = log10(wls)-log10(theMax/558);
	
	% Compute log optical density
	logDensity = a + ...
						   b*logWlsNorm.^2 + ...
							 c*logWlsNorm.^4 + ...
							 d*logWlsNorm.^6 + ...
							 e*logWlsNorm.^8 + ...
							 f*logWlsNorm.^10 + ...
							 g*logWlsNorm.^12 + ...
							 h*logWlsNorm.^14;
	T_absorbance(i,:) = 10.^logDensity;

    % If wavelength spacing is sufficiently fine, normalize explicitly
    % to a max of 1.  Seems like a good idea since this should always
    % be true of absorbance functions.  Don't want to do it if
    % wavelength spacing is coarse, because the peak may be between the
    % samples.
    if (S(2) <= 1)
        T_absorbance(i,:) = T_absorbance(i,:)/max(T_absorbance(i,:));
    end
end

