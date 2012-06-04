% RenderDemo
%
% Illustrates calibration interface for simple task of
% producing a uniform color patch of desired CIE xyY
% coordinates.
%
% The calculation is done with respect to the current
% PTB demonstration calibration file.
%
% The RGB values are gamma corrected and live in the range
% [0,1].  If they contain 0 or 1, the xyY coordinates
% requested may have been out of gamut.
%
% A uniform color patch is displayed in the MATLAB figure window.
% This is not a well-controlled display method, but does give
% a sense of the patch color if the calibration file is a
% reasonable description of the display.
%
% 4/26/97  dhb  Wrote it.
% 7/25/97  dhb  Better initialization.
% 3/12/98  dgp  Use Ask.
% 3/14/02  dhb  Update for OpenWindow.
% 4/03/02  awi  Merged in Windows changes.  On Windows we do not copy the result to the clipboard. 
% 4/13/02  awi	Changed "SetColorSpace" to new name "SetSensorColorSpace".
%				Changed "LinearToSettings" to new name "SensorToSettings".
% 12/21/02 dhb  Remove reliance on now obsolete OpenWindow/CloseWindow.
% 11/16/06 dhb  Start getting this to work with PTB-3.
% 11/22/06 dhb  Fixed except that Ask() needs to be fixed.
% 6/16/11  dhb  The PTB display section was out of date and didn't work.  I removed it.

% Clear out workspace
clear

% Load calibration file
cal = LoadCalFile('PTB3TestCal');
load T_xyz1931
T_xyz1931 = 683*T_xyz1931;
cal = SetSensorColorSpace(cal,T_xyz1931,S_xyz1931);
cal = SetGammaMethod(cal,0);

% Get xyY, render, and report.
xyY = input('Enter xyY (as a row vector) default [.3 .3 50]: ')';
if isempty(xyY)
    xyY = [.3 .3 50]';
end
XYZ = xyYToXYZ(xyY);
RGB = SensorToSettings(cal, XYZ);
fprintf('Computed RGB: [%g %g %g]\n', RGB(1), RGB(2), RGB(3));

% Make it an image
nX = 256; nY = 128;
theRGBCalFormat = RGB*ones(1,nX*nY);
theRGBImage = CalFormatToImage(theRGBCalFormat,nX,nY);

% Show in a Matlab figure window.  This will not be calibrated,
% but gives the general sense.  Use PTB display routines for
% more precise display.
figure; clf;
h = image(theRGBImage);
title('Here is the color');
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
set(gca,'XTick',[]);
set(gca,'YTick',[])

