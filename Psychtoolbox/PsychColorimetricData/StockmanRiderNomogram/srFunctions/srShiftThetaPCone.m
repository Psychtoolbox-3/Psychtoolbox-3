function [thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,wlShift,templatePeakWl)
% [thetaP,deltaThetaP] =  srShiftThetaPCone(thetaP,,wlShift,templatePeakWl)
%
% Shift thetaP values for a passed set of thetaP and the wlShift in nm, given template
% peak wavelength in nm.  This is implementing Equation 9 of the Stockman-Rider
% paper. This applies for photopigment thetaP's, not lens or macular pigment versions.

% History
%   2025-09-05  dhb  Modularize this calculation

% Equation 9 of the paper, once you work through the constant.
deltaThetaP = log10(templatePeakWl/(templatePeakWl+wlShift))/0.1187666467581842301;
thetaP= thetaP + deltaThetaP;

end