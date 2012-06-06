function blurshader = Create2DGaussianBlurShader(stddev, width)
% blurshader = Create2DGaussianBlurShader(stddev, width)
% EXPERIMENTAL - THIS FEATURE IS IN BETA-STAGE AND MAY CHANGE
% IN THE FUTURE OR NOT WORK AT ALL IN THE PRESENT!
%
% Create a GLSL shader which uses a 2D gaussian filter kernel to
% apply gaussian blur onto a texture. Returns the handle 'blurshader',
% which can be used via glUseProgram(blurshader) to activate the shader.
%
% stddev = Standard deviation of kernel to use. Defaults to 2.5.
% width  = Width of filter-kernel. Currently only a width of 5 is
% supported, so this argument is optional and defaults to 5.
%
% Currently you will need the Matlab imageprocessing toolbox for this to
% work, as we use the toolbox function 'fspecial' to compute the filter
% kernel.
%
% Example of use:
%
% Somewhere during initialization of your script you call:
% blurshader = Create2DGaussianBlurShader(2.5);
%
% When you want to draw an image blurred with a gaussian lowpass, you do:
%
% glUseProgram(blurshader);
% Screen('DrawTexture', win, mytexture, ....);
% glUseProgram(0);
%

% History:
% 15.04.2006 written by Mario Kleiner.

if nargin < 1
    stddev = 2.5;
end;

if nargin < 2 || isempty(width)
    width = 11;
end;

if width~=11
%    error('Create2DGaussianBlurShader: Invalid width specified. Currently only width=5 supported.');
end;

% Build and initialize gaussian blur shader (5x5, stddev=2.5):
try
    % Use the imageprocssing toolbox function fspecial for creation of
    % kernel:
    kernel = fspecial('gaussian', width, stddev);
catch
    error('Create2DGaussianBlurShader failed: You do not seem to have the image processing toolbox installed.');
end;

% Create shader from convolution kernel and return handle to it:
blurshader = EXPCreateStatic2DConvolutionShader(kernel);
return;
