function [minval,maxval] = MatMinMax(image)
% [min,max] = MatMinMax(image)
%
% Find the minimum and maximum values in a matrix.

% 12/27/04	dhb		Generalize to N-D input

minval = min(image(:));
maxval = max(image(:));
