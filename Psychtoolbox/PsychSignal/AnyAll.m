function [q] = AnyAll(in)
% [q] = AnyAll(in)
% returns true if any element of IN is non-zeros, false otherwise
% returns a scalar boolean. IN can be anything from a scalar to an N-D
% matrix

% DN 2008

q = any(in(:));
