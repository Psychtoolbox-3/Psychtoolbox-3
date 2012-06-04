function [vec,nRows,nCols] = ImageToVec(image)
% [vec,nRows,nCols] = ImageToVec(image)
%
% Take an image in matrix format and convert
% it to vector format.
%
% Also see VecToImage.
%
% 8/13/94		dhb		Added image size return


[nRows,nCols] = size(image);
vec = reshape(image,1,nRows*nCols)';
