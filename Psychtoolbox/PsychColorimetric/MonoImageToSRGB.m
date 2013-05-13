function rgbImage = MonoImageToSRGB(monoImage,xy,SCALE)
% rgbImage = MonoImageToSRGB(monoImage,xy,[SCALE])
%
% Take a single plane image in and convert to an sRGB standard
% color image out.  The chromaticity of the output image will
% be at the specify CIE x,y chromaticity value.
%
% Result is a uint8 image scaled into range [0,255].
%
% The value of SCALE is passed on to SRGBGammaCorrect, and
% there determines whether autoscaling is applied to the image.
%
% See MonoImageToSRGBTest.
%
% 6/15/11  dhb, ms  Wrote it.

% Set default
if (nargin < 3 || isempty(SCALE))
    SCALE = 1;
end

% First step.  Project input plane to CIE XYZ tristimulus image
% at specified chromaticity.
XYZ0 = xyYToXYZ([xy ; 1]);

% String out monochrome image values as a single row.
[monoImageCalFormat,nX,nY] = ImageToCalFormat(monoImage);

% Now convert to a 3 by n matrix, where each column is a scaled
% version of the target XYZ0 values
xyzImageUnscaledCalFormat = XYZ0*monoImageCalFormat;

% Convert to sRGB
rgbImagePrimaryCalFormat = XYZToSRGBPrimary(xyzImageUnscaledCalFormat);

% Gamma correct.  This scales the result so that the maximum of the
% RGB values comes out at 255.  Thus it is not so good for calibrated
% image processing, where we might want to preserve overall intensity
% differences.  But for just looking at stuff, it seems about right.
RGBImageCalFormat = uint8(SRGBGammaCorrect(rgbImagePrimaryCalFormat,SCALE));

% Turn it back into an image
rgbImage = CalFormatToImage(RGBImageCalFormat,nX,nY);

