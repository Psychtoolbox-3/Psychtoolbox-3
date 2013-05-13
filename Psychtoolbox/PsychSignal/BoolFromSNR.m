function q = BoolFromSNR(SNR,n)
% creates a boolean (selection) vector from an SNR (Signal-to-Noise Ratio)
% q = BoolFromSNR(SNR,n)

% 2008 DN  Wrote it.

% compute which proportion is signal for this SNR
prop    = (SNR)/(SNR+1);

% a percentage will be randomized to achieve the SNR
inds    = NRandPerm(n,round(n.*(1-prop)));
q       = ismember([1:n],inds);
