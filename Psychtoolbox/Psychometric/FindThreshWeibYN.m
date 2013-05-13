function threshold = FindThreshWeibYN(pyes, alpha, beta)
% threshold = FindThreshWeibYN(pcorrect, alpha, beta)
% 
%  Invert Weibull function to find threshold
%  for given pcorrect, alpha, and beta.
%  This function should invert ComputeWeibullYN().
%
%  threshold = alpha * pow( -1.0*log(1.0-pyes) , 1.0/beta );

  if (pyes>=0.0 && pyes<1.0)
    % ASSERTEND(pyes>=0.0 & pyes<1.0);
  end

  threshold = alpha .* ( -1.0 .* log(1.0-pyes)).^(1.0./beta);


