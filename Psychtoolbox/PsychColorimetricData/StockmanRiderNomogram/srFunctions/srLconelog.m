function y = srLconelog(nm, Lshift)
% y = srLconelog(nm, Lshift)
%
% This returns the fit to the Stockman-Sharp log10 L absorbance.  Paper
% Figure 1 and Table 1.
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog, srLserconelog

% History:
%   2025-09-05  dhb  Matlab first version as described above.

% Get Theta_P from wavelengths
thetaP = srWlsToThetaPCone(nm);

% Apply Lshift
Llmax_template = 551.9;
[thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,Lshift,Llmax_template);

% Table of coefficients.  Table 1 of paper, but more places.
% Extracted from the Python code by DHB because Claude
% didn't manage to do it.
c(1) = -42.9263580477;
c(2) = -2.0396798842;
c(3) = 75.9717833722;
c(4) = 57.3308210897;
c(5) = 6.5733913482;
c(6) = 8.1111028646;
c(7) = -38.7656494836;
c(8) = -21.4483453428;
c(9) = -5.9397465531;
c(10) = -3.3896198598;
c(11) = 9.5882997573;
c(12) = 3.2507563437;
c(13) = 1.4412770784;
c(14) = 0.3966003448;
c(15) = -0.7113921460;
c(16) = -0.0793542168;
c(17) = -0.0729797618;
c(18) = -0.0016552130;

% Compute the Fourier polynomial
y = c(1) + c(2)*cos(thetaP) + c(3)*sin(thetaP) + c(4)*cos(2*thetaP) + c(5)*sin(2*thetaP) + ...
    c(6)*cos(3*thetaP) + c(7)*sin(3*thetaP) + c(8)*cos(4*thetaP) + c(9)*sin(4*thetaP) + ...
    c(10)*cos(5*thetaP) + c(11)*sin(5*thetaP) + c(12)*cos(6*thetaP) + c(13)*sin(6*thetaP) + ...
    c(14)*cos(7*thetaP) + c(15)*sin(7*thetaP) + c(16)*cos(8*thetaP) + c(17)*sin(8*thetaP) + c(18);
end