function y = srSconelog(nm, Sshift)
% y = srMconelog(nm, Sshift)
%
% This returns the fit to the Stockman-Sharp log10 S absorbance.  Paper
% Figure 1 and Table 1.
%
% Adopted by Claude AI and DHB from Stockman-Rider paper and Python code.
%
% See also StockmanRiderDemo, srLMSconelog.

% History:
%   2025-09-05  dhb  Matlab first version as described above.

% Get Theta_P from wavelengths
thetaP = srWlsToThetaPCone(nm);

% Apply the wl shift as in Equation 1.
Slmax_template = 416.9;
[thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,Sshift,Slmax_template);

% Table of coefficients.  Table 1 of paper, but more places.
% Claude snagged these from the python code.
c = [207.3880950935, -6.3065623516, -393.7100478026, -315.6650602846, 19.2917535553, ...
    19.6414743488, 214.2211570447, 121.8584683485, -15.1820737886, -8.6774057156, ...
    -56.7596380441, -20.6318720369, 3.6934875040, 1.0483022480, 5.3656615075, ...
    0.7898783086, -0.1480357836, 0.0002358232];

% Compute the Fourier polynomial
y = c(1) + c(2)*cos(thetaP) + c(3)*sin(thetaP) + c(4)*cos(2*thetaP) + c(5)*sin(2*thetaP) + ...
    c(6)*cos(3*thetaP) + c(7)*sin(3*thetaP) + c(8)*cos(4*thetaP) + c(9)*sin(4*thetaP) + ...
    c(10)*cos(5*thetaP) + c(11)*sin(5*thetaP) + c(12)*cos(6*thetaP) + c(13)*sin(6*thetaP) + ...
    c(14)*cos(7*thetaP) + c(15)*sin(7*thetaP) + c(16)*cos(8*thetaP) + c(17)*sin(8*thetaP) + c(18);
end