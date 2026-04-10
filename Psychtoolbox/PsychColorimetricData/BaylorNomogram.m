function T_absorbance = BaylorNomogram(S,lambdaMax)
% T_absorbance = BaylorNomogram(S,lambdaMax)
%
% Compute spectral sensitivities according to the
% nomogram provided in Baylor, Nunn, and Schnapf, 1987.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to input incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 6/22/96  dhb  Wrote it.
% 10/16/97 dhb  Add comment about energy units.

% These are the coefficients for the polynomial 
% approximation.
aN = [-5.2734 -87.403 1228.4 -3346.3 -5070.3 30881 -31607];

% Get wls argument.
S = MakeItS(S);
wls = MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absorbance = zeros(nT,nWls);
wlsum = wls/1000;

for i = 1:nT
    wlsVec = log10( (1 ./ wlsum)*lambdaMax(i)/561)';
    logS = aN(1) + aN(2)*wlsVec + aN(3)*wlsVec.^2 + aN(4)*wlsVec.^3 + ...
        aN(5)*wlsVec.^4 + aN(6)*wlsVec.^5 + aN(7)*wlsVec.^6;
    T_absorbance(i,:) = 10.^logS;

    % If wavelength spacing is sufficiently fine, normalize explicitly
    % to a max of 1.  Seems like a good idea since this should always
    % be true of absorbance functions.  Don't want to do it if
    % wavelength spacing is coarse, because the peak may be between the
    % samples.
    if (S(2) <= 1)
        T_absorbance(i,:) = T_absorbance(i,:)/max(T_absorbance(i,:));
    end
end
