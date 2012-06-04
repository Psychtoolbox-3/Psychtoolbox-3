function trolands = PowerToTrolands(wl,power,S_vLambda,vLambda)
% trolands = PowerToTrolands(wl,power,S_vLambda,vLambda)
%
% Convert from power (in watts/deg-2) to trolands for a
% monochromatic light.  Wavelength specified in nm,
% vLambda should be CIE photopic luminosity function.
%
% See W&S, p. 103.
%
% 8/15/96  dhb, abp  Wrote the header.

wls_vLambda = MakeItWls(S_vLambda);
index = find(wl == wls_vLambda);
trolands = power*vLambda(index)*2.242e12;
