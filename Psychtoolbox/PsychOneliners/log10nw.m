function y=log10nw(x)
% y=log10nw(x)
% "nw" stands for "No Warnings". Equivalent to Matlab's built-in LOG10 function, but suppresses warnings, 
% e.g. "Warning: Log of zero."
%
%   LOG10NW(X) is the base 10 logarithm of the elements of X.   
%   Complex results are produced if X is negative.
if nargin~=1
	error('Usage: y=log10nw(x)')
end
oldWarning=warning;
warning('off');
y = log(x) ./ log(10);
warning(oldWarning);
	
