function newSourceMat=AlphaSourceTerm(sourceFactorStr, sourceMat, destinationMat)

% newSourceMat=AlphaSourceTerm(sourceFactorStr, sourceMat, destinationMat)
%
% AlphaSourceTerm simuluates the step of multplying an alpha factor with the
% source image in OpenGL alpha blending.    
%
% Apply an OpenGL source factor rule such as 'GL_ZERO' to the image matrix
% 'sourceMat', returing the dot product of sourceMat and the alpha factor
% selected by string 'sourceFactorStr'.  
%
% The destination image matrix 'destinationMat' is required because source
% factor strings may select alpha values from the destination.
%
% AlphaSourceTerm calculates with double-precision (64-bit) floating point
% arithmatitic whereas the precision of an OpenGL renderer which it
% simulates is unspecified, except that OpenGL guarantees perfect precision
% for alpha values 0 and 1 (255 via Screen).  Comparison of AlphaSourceTerm
% with the OpenGL renderer shows that alpha multiplicaion discards up to
% one bit of precision.  
%
% see also: AlphaDestinationProduct PsychAlphaBlending

% HISTORY
% 
% mm/dd/yy
% 
%  2/11/05  awi wrote it.


sourceMatDims=size(sourceMat);
destinationMatDims=size(destinationMat);
if length(sourceMatDims) ~= 3 ||  sourceMatDims(3) ~= 4                     
    error('Argument ''sourceMat'' must be a matrix with dimensions [x,y,4]');
end
if length(destinationMatDims) ~= 3 ||  destinationMatDims(3) ~= 4                     
    error('Argument ''destinationMat'' must be a matrix with dimensions [x,y,4]');
end
sourceMatAlphaPlane=sourceMat(:,:,4);
destinationMatAlphaPlane=destinationMat(:,:,4);
    

if strcmp(sourceFactorStr, 'GL_ZERO')
    newSourceMat=AlphaTimes(sourceMat, 0);
elseif strcmp(sourceFactorStr, 'GL_ONE')
    newSourceMat=AlphaTimes(sourceMat, 255);
elseif strcmp(sourceFactorStr, 'GL_DST_COLOR')
    newSourceMat=AlphaTimes(sourceMat, destinationMat);
elseif strcmp(sourceFactorStr, 'GL_ONE_MINUS_DST_COLOR')
    newSourceMat=AlphaTimes(sourceMat, 255 - destinationMat);
elseif strcmp(sourceFactorStr, 'GL_SRC_ALPHA');
    newSourceMat=AlphaTimes(sourceMat, sourceMat(:,:,4));
elseif strcmp(sourceFactorStr, 'GL_ONE_MINUS_SRC_ALPHA')
    newSourceMat=AlphaTimes(sourceMat, 255 - sourceMat(:,:,4));
elseif strcmp(sourceFactorStr, 'GL_DST_ALPHA');
    newSourceMat=AlphaTimes(sourceMat, destinationMat(:,:,4));
elseif strcmp(sourceFactorStr, 'GL_ONE_MINUS_DST_ALPHA');
    newSourceMat=AlphaTimes(sourceMat, 255-destinationMat(:,:,4));
elseif strcmp(sourceFactorStr, 'GL_SRC_ALPHA_SATURATE');
    minAlphaPlane=min(sourceMat(:,:,4), 255-destinationMat(:,:,4));
    newSourceMat=AlphaTimes(sourceMat, minAlphaPlane);
    newSourceMat(:,:,4)=ones(sourceMatDims(1), sourceMatDims(2)) * 255; 
else
    error('Argument "sourceFactorStr" is unrecognized or invalid');
end
    

