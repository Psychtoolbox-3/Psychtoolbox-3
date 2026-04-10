function T_absorbance = LambNomogram(S,lambdaMax)
% T_absorbance = LambNomogram(S,lambdaMax)
%
% Compute spectral sensitivities according to the
% nomogram provided in Lamb, 1995, Vision Research,
% Vol. 35, pp. 3083-3091, equation 2'.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 8/20/98 dhb  Wrote it.

% These are the coefficients for Equation 2.
a = 70; b = 28.5; c = -14.1;
A = 0.880; B = 0.924; C = 1.104; D = 0.655;

% Get wls argument.
S = MakeItS(S);
wls = MakeItWls(S);
[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absorbance = zeros(nT,nWls);

for i = 1:nT
	theMax = lambdaMax(i);
	wlarg = theMax ./ wls';
	T_absorbance(i,:) = 1 ./ ( ...
					  exp( a*(A-wlarg) ) + ...
						exp( b*(B-wlarg) ) + ...
						exp( c*(C-wlarg) ) + ...
						D ...
						);

    % If wavelength spacing is sufficiently fine, normalize explicitly
    % to a max of 1.  Seems like a good idea since this should always
    % be true of absorbance functions.  Don't want to do it if
    % wavelength spacing is coarse, because the peak may be between the
    % samples.
    if (S(2) <= 1)
        T_absorbance(i,:) = T_absorbance(i,:)/max(T_absorbance(i,:));
    end
end
