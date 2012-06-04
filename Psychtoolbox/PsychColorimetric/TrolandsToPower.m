function power = TrolandsToPower(wl,trolands,S_vLambda,vLambda)
% power = TrolandsToPower(wl,trolands,S_vLambda,vLambda)
%
% Convert from photopic trolands to power (watts/degree-2)
% for a monochromatic light.
%
% See W&S, p. 103.
%
% 8/15/96  dhb, abp  Wrote the header.

wls_vLambda = MakeItWls(S_vLambda);
index = find(wl == wls_vLambda);
power = trolands/(vLambda(index)*2.242e12);
