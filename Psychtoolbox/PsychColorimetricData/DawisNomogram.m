function [T_absorbance] = DawisNomogram(S,lambdaMax)
% [T_absorbance] = DawisNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided in Dawis, 1981, Vision Research,
% Vol. 21, pp. 1427-1430.
%
% T_absorbance contains the absorbance.
% This is specified according to the CVRL web page
% definition absorbance = log(I_incident/I_transmitted),
% so the numbers are all positive.  The peak is normalized
% to one.
%
% For low density, the absorbance has the same shape as
% the spectral senstivity.
%
% The nomogram is shifted along the wavelength axis
% using a multiplicative rather than additive procedure.
%
% The result is in quantal units, in the sense that to compute
% absorbtions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 10/30/97 dhb  Wrote it.
% 07/01/03 dhb  Add computation of T_absorbance.

% These are the coefficients for the polynomial 
% approximation, taken from Table 1.  We implement
% the A1-based pigment nomogram, not the A2.
Lmax = [432 ; 502 ; 562];
bN = [0.346325  -35.1001  63.3807 -125.466  10962.7  -16244.8 -210671.0   -23776.9 ; ...
     -0.0106836 -28.28   148.133  -498.627  -1457.94  12799.4    -789.371 -60749.2 ; ...
		 -0.228262  -22.9974  87.0027 -636.336   2624.45   4948.6  -45944.6    65688.5 ];
lmaxLow =  [410 ; 470 ; 530];
lmaxHigh = [470 ; 530 ; 610];
L1 = [380 ; 400 ; 430 ];
L2 = [510 ; 620 ; 690 ];

% Get wls argument.
wls = MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absorbance = zeros(nT,nWls);
wlsum = wls/1000;

for i = 1:nT
	theMax = lambdaMax(i);
	if (theMax >= lmaxLow(1) && theMax <= lmaxHigh(1))
		which = 1;
	elseif (theMax > lmaxLow(2) && theMax <= lmaxHigh(2))
		which = 2;
	elseif (theMax > lmaxLow(3) && theMax <= lmaxHigh(3))
		which = 3;
	else
		error(sprintf('Lambda Max %g not in range of nomogram\n',theMax));
	end
	wlsVec = (theMax ./ wls') - 1;
	logS = zeros(1,nWls);
	for k = 1:8
		logS = logS + bN(which,k)*wlsVec.^k;
	end
	T_absorbance(i,:) = logS;
	
	% Zero sensitivity outsize valid range.  I shift the
	% range in the table to slide multiplicatively with 
	% theMax.
	zeroLow = theMax*L1(which)/Lmax(which);
	zeroHigh = theMax*L2(which)/Lmax(which);
	index = find(wls < zeroLow);
	T_absorbance(i,index) = -Inf;
	index = find(wls > zeroHigh);
	T_absorbance(i,index) = -Inf;
end
T_absorbance = 10.^T_absorbance;
