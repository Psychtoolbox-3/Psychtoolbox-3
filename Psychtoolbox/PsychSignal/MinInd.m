function [m,inds] = MinInd(mat)
% [m,inds] = MinInd(mat)
% returns minimum of an N-D matrix and the indices at which this minimum
% occurs
% indices are returned in a [#mins x ndims] matrix

% DN 2008-07-21 wrote it

m                   = min(mat(:));
inds                = FindInd(mat==m);
