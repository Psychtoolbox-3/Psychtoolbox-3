function [red,green,blue] = UnpackColorImage(image)	
% [red,green,blue] = UnpackColorImage(image)	
% 
% Take a color image and unpack it into three
% three image planes.
% 
% Particularly useful for fixing old code
% that used SCREEN('GetColorImage',....)
%
% 11/24/02  jmh, dhb  Wrote it.

red = image(:,:,1);
green = image(:,:,2);
blue = image(:,:,3);

