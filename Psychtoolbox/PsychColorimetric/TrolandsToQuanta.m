function quanta = TrolandsToQuanta(wl,trolands,S_vLambda,vLambda)
% quanta = TrolandsToQuanta(wl,trolands,S_vLambda,vLambda)
%
% Convert from photopic trolands to quanta (sec-1 degree-2)
% for a monochromatic light.
%
% See W&S, p. 103.
%
% 8/15/96  dhb, abp  Wrote the header.
% 8/22/96  dhb       Wavelength should be in meters.

wls_vLambda = MakeItWls(S_vLambda);
index = find(wl == wls_vLambda);
quanta = trolands*(wl*1e-9)/(vLambda(index)*4.454e-13);
