function [maxval] = MatMax(image)
% [maxval] = MatMax(image)
% Find the maximum value in a matrix.

% 6/13/12    dn     Generalize to N-D input by forwarding to MaxAll()

maxval = MaxAll(image);
