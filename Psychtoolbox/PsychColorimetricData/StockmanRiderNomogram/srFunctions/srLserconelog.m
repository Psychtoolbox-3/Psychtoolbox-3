function y = srLserconelog(nm, Lshift)
% y = srLserconelog(nm, Lshift)
%
% Returns the codon 180 serine variant L cone photopigment absorbance Fourier fit.
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog, srLconelog

% History:
%   2025-09-05  dhb  Matlab first version as described above.

% Get Theta_P from wavelengths
thetaP = srWlsToThetaPCone(nm);

% Apply Lshift
Lserlmax_template = 553.1;
[thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,Lshift,Lserlmax_template);

% Table of coefficients.  Extracted from Python code by Claude AI.  In Table 4
% of the paper.
c = [-42.417608560, -2.656791612, 75.011093607, 56.477062776, 7.509397607, ...
    9.061442173, -38.068488495, -20.974610259, -6.642746250, -3.785039126, ...
    9.322071459, 3.134494745, 1.603799055, 0.439302358, -0.676958684, ...
    -0.072988371, -0.078857510, -0.004264105];

% Bada-bang, bada-bing
y = c(1) + c(2)*cos(thetaP) + c(3)*sin(thetaP) + c(4)*cos(2*thetaP) + c(5)*sin(2*thetaP) + ...
    c(6)*cos(3*thetaP) + c(7)*sin(3*thetaP) + c(8)*cos(4*thetaP) + c(9)*sin(4*thetaP) + ...
    c(10)*cos(5*thetaP) + c(11)*sin(5*thetaP) + c(12)*cos(6*thetaP) + c(13)*sin(6*thetaP) + ...
    c(14)*cos(7*thetaP) + c(15)*sin(7*thetaP) + c(16)*cos(8*thetaP) + c(17)*sin(8*thetaP) + c(18);

end