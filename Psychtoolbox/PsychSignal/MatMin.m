function [minval] = MatMin(image)
% [minval] = MatMin(image)
% Find the minumum value in a matrix.

% 6/13/12    dn     Generalize to N-D input by forwarding to MinAll()

minval = MinAll(image);
