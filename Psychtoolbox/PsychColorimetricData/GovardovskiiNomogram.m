function T_absorbance = GovardovskiiNomogram(S,lambdaMax)
% T_absorbance = GovardovskiiNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided in Victor I. Govardovskii et al., 2000,
% Visual Neuroscience, Vol. 17, pp. 509-528.
%
% The polynomial approximation has two bands.
% Alpha-band, the major band of A1 pigments is taken from
% Equation (1) and equation (2) in wavelength range [350,700] nm
% (see page 515 in paper)
%
% Beta-band, the minor band of A1 pigments is taken from
% Equation (3) and equation (5a, 5b) in wavelgnth range [350,700] nm
% (see page 516 in paper)
%
% T_absorbance contains the absorbance.
% This is specified according to the CVRL web page
% definition absorbance = log(I_incident/I_transmitted),
% so the numbers are all positive.  The peak is normalized
% to one.
%
% For low density, the absorbance has the same shape as
% the spectral senstivity.
%
% The result is in quantal units, in the sense that to compute
% absorbtions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 03/08/2002 ly  Wrote it starting from DawisNomogram.

% Valid range of wavelengh for A1-based visual pigments
% (see page 516 in paper)
Lmin = 330;
Lmax = 700;

% Valid range of lambdamax value
% (taken from figure 4, page 514 in paper)
lmaxLow = 350;
lmaxHigh = 600;

% Alpha-band parameters
% A = 69.7, B = 28, C = 14.9, D = 0.674 for equation (1) (page 515)
A_B_C = [69.7, 28, -14.9];
D = 0.674;

% b = 0.922, c = 1.104 for equation (1) (page 515)
b_c = [0.922, 1.104];

% Beta-band parameters
% Abeta = 0.26 for equation (4) (page 516)
Abeta = 0.26;

% Get wls argument.
wls = MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absorbance = zeros(nT,nWls);

for i = 1:nT
    theMax = lambdaMax(i);
    if (theMax > lmaxLow && theMax < lmaxHigh)
        
        % alpha-band polynomial
        %
        % Parameter a depends on lambdamax, see equation (2) (page 515)
        a = 0.8795 + 0.0459*exp(-(theMax-300)^2/11940);
        a_b_c = [a, b_c];
        
        x = theMax./wls;
        
        % midStepN, N = 1, 2, ... are the middle steps in the caculation.
        midStep1 = exp (ones(nWls,1)*(A_B_C.*a_b_c) - x*A_B_C);
        midStep2 = sum(midStep1,2) + D;
        
        % Result of equation (1) (page 515)
        S_x = 1./midStep2;
        
        % Beta-band polynomial
        
        % Parameter bbeta depends on lambdamax, see equation (5b) (page 516)
        bbeta = -40.5 + 0.195*theMax;
        
        % Conversion of lambdamax to parameter lambdaMaxbeta, see equation (5a) (page 516)
        lambdaMaxbeta = 189 + 0.315*theMax;
        
        % midStepN, N = 1, 2, ... are the middle steps in the caculation.
        midStep1 = -((wls - lambdaMaxbeta * ones (nWls,1)) / bbeta).^2;
        midStep2 = Abeta * exp (midStep1);
        
        % Result of equation (4) (page 516)
        S_beta = midStep2;
        
        % alpha band and beta band together.
        T_absorbance(i,:) = (S_x + S_beta)';
        
        % Zero sensitivity outsize valid range.
        index = find(wls < Lmin);
        T_absorbance(i,index) = zeros(size(index))';
        index = find(wls > Lmax);
        T_absorbance(i,index) = zeros(size(index))';
        
    else
        error(sprintf('Lambda Max %g not in range of nomogram\n',theMax));
    end
    
end
