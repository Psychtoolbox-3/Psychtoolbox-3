function trolands = QuantaToTrolands(wl,quanta,S_vLambda,vLambda)
% trolands = QuantaToTrolands(wl,quanta,S_vLambda,vLambda)
%
% Convert from quanta (sec-1 deg-2) to to trolands for a
% monochromatic light.
%
% See W&S, p. 103.
%
% 8/15/96  dhb, abp  Wrote the header.
% 8/22/96  dhb       Wavelength should be in meters.

wls_vLambda = MakeItWls(S_vLambda);
index = find(wl == wls_vLambda);
trolands = (quanta/(wl*1e-9))*vLambda(index)*4.454e-13;
