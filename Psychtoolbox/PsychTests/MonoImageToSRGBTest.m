% MonoImageToSRGBTest
%
% Test that our litte routine that converts a monochrome image to a
% color image works about right.
%
% 6/15/11  dhb, ms  Wrote it.

%% Clear
clear; close all;

%% Set chromaticity.  This
% basically sets the "tint' of
% the produced image.
xy = [0.3 0.35]';

%% Read in a monochrome test image
testImageName = 'konijntjes1024x768gray.jpg';
if (~exist(testImageName,'file'))
   error('Someone removed the test image from the Psychtoolbox.  Find a new one.');
end
testImage = double(imread(testImageName));


%% Call our routine
rgbImage = MonoImageToSRGB(testImage,xy);

%% Display
figure;
subplot(1,2,1);
imshow(testImage/max(testImage(:)));
subplot(1,2,2);
imshow(rgbImage);
