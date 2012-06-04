function image = VecToImage(vec,nRows,nCols)
% image = VecToImage(vec,[nRows,nCols])
%
% Convert an image from vector to matrix format.
% Image must be square for this to work properly.
%
% Also see ImageToVec.
%
% 8/13/94		dhb		Added optional nRows,nCols arguments.

% Figure out image size
[m,n] = size(vec);

if (nargin == 1)
	nRows = floor(sqrt(m));
	nCols = nRows;
elseif (nRows*nCols ~= m)
	error('Passed image dimensions do not match vector size');
end

% Reshape it into an image
image = reshape(vec,nRows,nCols);

