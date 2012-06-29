function [vec,nRows,nCols] = ImageToVec(image)
% [vec,nRows,nCols] = ImageToVec(image)
%
% Take an image in matrix format and convert
% it to vector format.
%
% Also see VecToImage.
%
% 8/13/94		dhb		Added image size return
% 6/13/12		 dn		No need for reshape call


[nRows,nCols] = size(image);
vec = image(:);
