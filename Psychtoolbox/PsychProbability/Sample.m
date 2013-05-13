function [Y,index] = Sample(X)
% [Y,index] = Sample(X)
%
% Returns a random sample from X.
% If X is a vector, returns a scalar sample from X, so Y = X(index).
% If X is an m-by-n matrix, m>1 && n>1, samples each column of X:
%	for j=1:n, Y(j)=X(index(j),j).
%
% Also see RandSample, Shuffle, SORT, and Randi.

% 6/29/96	dgp	    Wrote it.
% 7/26/04   awi     Changed 'randi' to 'Randi'.  Matlab 7 is case-sensitive.    
% 1/29/05	dgp	    Cosmetic.

[n,m]=size(X);
if n==1 || m==1
	index=Randi(length(X));
	Y=X(index);
else
	index=Randi(n*ones(1,m));
	Y=X(index+(0:m-1)*n);
end
