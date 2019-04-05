function T = BaylorNomogram(S,lambdaMax)
% T = BaylorNomogram(S,lambdaMax)
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
wls = MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T = zeros(nT,nWls);
wlsum = wls/1000;

for i = 1:nT
	wlsVec = log10( (1 ./ wlsum)*lambdaMax(i)/561)';
	logS = aN(1) + aN(2)*wlsVec + aN(3)*wlsVec.^2 + aN(4)*wlsVec.^3 + ...
					 aN(5)*wlsVec.^4 + aN(6)*wlsVec.^5 + aN(7)*wlsVec.^6;
	T(i,:) = 10.^logS;
end
