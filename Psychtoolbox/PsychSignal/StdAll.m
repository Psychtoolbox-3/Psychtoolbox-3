function s = StdAll(in)
% [s] = StdAll(in)
% returns the standard deviation of input IN. IN can be anything from a
% scalar to an N-D matrix
%
% DN 2008

s = std(in(:));
