function out = AltRound(nr,factor)
% out = AltRound(nr,factor)
% rounds the elements in NR to the nearest multiple
% of FACTOR.
% 
% FACTOR must be either a scalar or of the same shape as NR.
%
% Thus:
%   AltRound(12.26,.25) = 12.25;
%   AltRound(12.2 ,0.6) = 12,
%   but also
%   AltRound(1449 ,100) = 1400.
%
% DN 2008
% DN 2009-02-02 Support non-scalar factor input

psychassert(all(isnumeric(factor(:))) && all(isfinite(factor(:))),'Factor argument must be a numeric and finite.')
psychassert(isscalar(factor) || isequal(size(factor),size(nr)),'Factor argument must be scalar or of the same shape as nr argument.')


factor = 1./factor;

out = round(nr.*factor)./factor;
