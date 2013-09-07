function threshold = FindThreshWeibTAFC(pCorrect,alpha,beta)
% threshold = FindThreshWeibTAFC(pCorrect,alpha,beta)
%
% Invert Weibull function to find threshold for given pcorrect,
% alpha, and beta.
%
% This function should invert ComputeWeibTAFC().
%
% 8/25/94 dhb, ccc    New name, fixed comments

if (pCorrect <= 0.0 || pCorrect > 1.0)
  error('pCorrect is out of bounds');
end

threshold = alpha .* ( -1.0*log(2.0 .* (1.0-pCorrect))) .^(1.0./beta);

