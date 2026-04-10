function thetaP = srWlsToThetaPCone(nm)
% thetaP = srWlsToThetaPCone(nm)(nm)
%
% Get thetap from wavelength in nm, according to Eq 2 of the Stockman-Rider paper.  This
% applies to photopigment thetaPs, not to lens and macular pigment versions.
%
% See also StockmanRiderDemo, srLconelog, srMconelog, srSconelog

% History
%   2025-09-05  dhb  Modularized this little bit.

% This is Eq 2 of the paper, with x representing theta_p.
%   log10(360) = 2.5563
%   1/(pi/log10(850/360) = 0.1188
% This puts in more places than display in Matlab.
thetaP = (log10(nm)-2.556302500767287267)/0.1187666467581842301;

end