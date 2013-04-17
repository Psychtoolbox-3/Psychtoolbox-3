function luminance = TrolandsToLum(trolands,pupilAreaMM)
% luminance = TrolandsToLum(trolands,pupilAreaMM)
%
% Convert photopic/scotopic trolands to corresponding cd/m2.
%
% 7/29/03  dhb  Wrote it.

luminance = trolands/pupilAreaMM;
