function [Y,index] = Shuffle(X)
% [Y,index] = Shuffle(X)
%
% Randomly sorts X.
% If X is a vector, sorts all of X, so Y = X(index).
% If X is an m-by-n matrix, sorts each column of X, so
%	for j=1:n, Y(:,j)=X(index(:,j),j).
%
% Also see SORT, Sample, Randi, and RandSample.

% xx/xx/92  dhb  Wrote it.
% 10/25/93  dhb  Return index.
% 5/25/96   dgp  Made consistent with sort and "for i=Shuffle(1:10)"
% 6/29/96	  dgp  Edited comments above.
% 5/18/02   dhb  Modified code to do what comments say, for matrices.
% 6/2/02    dhb  Fixed bug introduced 5/18.

[null,index] = sort(rand(size(X)));
[n,m] = size(X);
Y = zeros(size(X));
if n == 1 || m == 1
	Y = X(index);
else
	for j = 1:m
		Y(:,j)  = X(index(:,j),j);
	end
end
 
