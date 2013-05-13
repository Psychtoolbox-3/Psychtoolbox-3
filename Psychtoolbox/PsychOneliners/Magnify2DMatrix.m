function destination = Magnify2DMatrix(source, scalingFactor) 
% destination = Magnify2DMatrix(sourceMatrix, scalingFactor)
%
% Magnifies a 2-D, 3-D or N-D 'sourceMatrix' by a factor specified by
% 'scalingFactor'. Size of 3rd and higher dimensions will not be scaled.
%

% 10/15/06 rhh Wrote it using lots of loops.
% 11/12/06 rhh Revised it.  Denis Pelli suggested revising the code to take
%               advantage of Matlab's matrix processing abilities and David Brainard
%               showed explicitly how this could be done.
% 05/09/08 DN  generated copy instruction indices with a different method,
%               speeding up this function. Added support for 3D matrices
% 13/06/12 DN  Now support input with arbitrary number of dimensions. Moved
%               implementation to Expand()

destination = Expand(source, scalingFactor);
