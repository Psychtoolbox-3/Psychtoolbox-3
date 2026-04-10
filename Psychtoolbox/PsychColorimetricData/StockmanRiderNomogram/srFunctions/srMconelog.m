function y = srMconelog(nm, Mshift)
% y = srMconelog(nm, Mshift)
%
% This returns the fit to the Stockman-Sharp log10 M absorbance.  Paper
% Figure 1 and Table 1.
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog.

% History:
%   2025-09-05  dhb  Matlab first version as described above.

% Get Theta_P from wavelengths
thetaP = srWlsToThetaPCone(nm);

% Apply the wl shift as in Equation.  Mlmax given as 529.8 in the paper,
% and that matches the polynomial when it is computed for 0 shift at 0.1
% nm spacing.  This value is 529.9 in the python code, but I'm going with
% the paper.  0.1 nm is not of any practical significance, and this only 
% has a second order effect when you shift because it determines how
% the thetaPs get shifted, and it won't be much different.
Mlmax_template = 529.8;
[thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,Mshift,Mlmax_template);

% Table of coefficients.  Table 1 of paper, but more places.
% Claude snagged these from the python code.
c = [-210.6568853069, -0.1458073553, 386.7319763250, 305.4710584670, 5.0218382813, ...
    6.8386224350, -208.2062335724, -118.4890200521, -5.7625866330, -3.7973553168, ...
    55.1803460639, 19.9728512548, 1.8990456325, 0.6913410864, -5.0891806213, ...
    -0.7070689492, -0.1419926703, 0.0005894876];

% Compute the Fourier polynomial
y = c(1) + c(2)*cos(thetaP) + c(3)*sin(thetaP) + c(4)*cos(2*thetaP) + c(5)*sin(2*thetaP) + ...
    c(6)*cos(3*thetaP) + c(7)*sin(3*thetaP) + c(8)*cos(4*thetaP) + c(9)*sin(4*thetaP) + ...
    c(10)*cos(5*thetaP) + c(11)*sin(5*thetaP) + c(12)*cos(6*thetaP) + c(13)*sin(6*thetaP) + ...
    c(14)*cos(7*thetaP) + c(15)*sin(7*thetaP) + c(16)*cos(8*thetaP) + c(17)*sin(8*thetaP) + c(18);
end