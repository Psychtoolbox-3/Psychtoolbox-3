function alpha = FindAlphaWeibTAFC(pCorrect,threshold,beta)
% threshold = FindThreshWeibTAFC(pCorrect,threshold,beta)
%
% Find alpha so that the Weibull function with passed beta
% predicts the pCorrect threshold to be the passed value
% of threshold.
%
% This computation is mostly useful for testing.
%
% 8/25/94		dhb, ccc		Wrote it.

if (pCorrect <= 0.0 || pCorrect > 1.0)
  error('pCorrect is out of bounds');
end

% Find the alpha.  We didn't really think about this, we just
% assumed that the code in FindThreshWeibTAFC was correct and
% went from there.
alpha = threshold ./ ( -1.0*log(2.0 .* (1.0-pCorrect))) .^(1.0./beta);

