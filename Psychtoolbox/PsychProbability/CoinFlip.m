function x = CoinFlip(ndraws,p)
% x = CoinFlip(ndraws,p)
%
% Generate a list x of zeros and ones according to coin flipping (i.e.
% Bernoulli) statistics with probability p of getting a 1.

% 1/20/97  dhb  Delete obsolet rand('uniform').
% 7/24/04  awi  Cosmetic.
% 6/13/12   dn  Use boolean vectors, simplifies code


% Generate ndraws random variables on the real
% interval [0,1).
unif = rand(ndraws,1);

% Find all of the ones that are less than p.
% On average, this proportion will be p.
index = unif < p;

% Generate an array of zeros and then set
% the ones found in the previous step to 1.
x = zeros(ndraws,1);
x(index) = 1;



