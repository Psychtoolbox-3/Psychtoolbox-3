% DKLDemo
%
% Demonstrates the DKL color space routines.
% Produces a picture of the DKL isoluminant plane
% in a figure window.  This picture is not calibrated.
%
% You can chooose the image size and what cone/luminance
% data you want to use to define the space.
%
% 4/9/05	dhb		Wrote it.
% 5/26/11   dhb     Updated to use PTB3TestCal, SetSensorColorSpace
%           dhb     Remove OS9 option.

% Clear
clear; close all

% Define parameters.  Image size should be an odd number.
imageSize = 513;
whichCones = 'SmithPokorny';

% Load spectral data and set calibration file
switch (whichCones)
	case 'SmithPokorny',
		load T_cones_sp
		load T_xyzJuddVos
		S_cones = S_cones_sp;
		T_cones = T_cones_sp;
		T_Y = 683*T_xyzJuddVos(2,:);
		S_Y = S_xyzJuddVos;
		T_Y = SplineCmf(S_Y,T_Y,S_cones);
	case 'StockmanSharpe'
		load T_cones_ss2
		load T_ss2000_Y2
		S_cones = S_cones_ss2;
		T_cones = T_cones_ss2;
		T_Y = 683*T_ss2000_Y2;
		S_Y = S_ss2000_Y2;
		T_Y = SplineCmf(S_Y,T_Y,S_cones);
end
cal = LoadCalFile('PTB3TestCal');
calLMS = SetSensorColorSpace(cal,T_cones,S_cones);
calLMS = SetGammaMethod(calLMS,1);
calLum = SetSensorColorSpace(cal,T_Y,S_Y);

% Define background.  Here we just take the
% monitor mid-point.
bgLMS = PrimaryToSensor(calLMS,[0.5 0.5 0.5]');

% Get matrix that transforms between incremental
% cone coordinates and DKL coordinates 
% (Lum, RG, S).
M_ConeIncToDKL = ComputeDKL_M(bgLMS,T_cones,T_Y);
M_DKLToConeInc = inv(M_ConeIncToDKL);

% Now find incremental cone directions corresponding
% to DKL isoluminant directions.
rgConeInc = M_DKLToConeInc*[0 1 0]';
sConeInc = M_DKLToConeInc*[0 0 1]';

% These directions are not scaled in an interesting way,
% need to scale them.  Here we'll find units so that 
% a unit excursion in the two directions brings us to
% the edge of the monitor gamut, with a little headroom.
bgPrimary = SensorToPrimary(calLMS,bgLMS);
rgPrimaryInc = SensorToPrimary(calLMS,rgConeInc+bgLMS)-bgPrimary;
sPrimaryInc = SensorToPrimary(calLMS,sConeInc+bgLMS)-bgPrimary;
rgScale = MaximizeGamutContrast(rgPrimaryInc,bgPrimary);
sScale = MaximizeGamutContrast(sPrimaryInc,bgPrimary);
rgConeInc = 0.95*rgScale*rgConeInc;
sConeInc = 0.95*sScale*sConeInc;

% If we find the RGB values corresponding to unit excursions
% in rg and s directions, we should find a) that the luminance
% of each is the same and b) that they are all within gamut.
% In gamut means that the primary coordinates are all bounded
% within [0,1].
rgPlusLMS = bgLMS+rgConeInc;
rgMinusLMS = bgLMS-rgConeInc;
sPlusLMS = bgLMS+sConeInc;
sMinusLMS = bgLMS-sConeInc;
rgPlusPrimary = SensorToPrimary(calLMS,rgPlusLMS);
rgMinusPrimary = SensorToPrimary(calLMS,rgMinusLMS);
sPlusPrimary = SensorToPrimary(calLMS,sPlusLMS);
sMinusPrimary = SensorToPrimary(calLMS,sMinusLMS);
if (any([rgPlusPrimary(:) ; rgMinusPrimary(:) ; ...
		sPlusPrimary(:) ; sMinusPrimary(:)] < 0))
	fprintf('Something out of gamut low that shouldn''t be.\n');
end
if (any([rgPlusPrimary(:) ; rgMinusPrimary(:) ; ...
		sPlusPrimary(:) ; sMinusPrimary(:)] > 1))
	fprintf('Something out of gamut high that shouldn''t be.\n');
end
bgLum = PrimaryToSensor(calLum,bgPrimary);
rgPlusLum = PrimaryToSensor(calLum,rgPlusPrimary);
rgMinusLum = PrimaryToSensor(calLum,rgMinusPrimary);
sPlusLum = PrimaryToSensor(calLum,sPlusPrimary);
sMinusLum = PrimaryToSensor(calLum,sMinusPrimary);
lums = sort([bgLum rgPlusLum rgMinusLum sPlusLum sMinusLum]);
fprintf('Luminance range in isoluminant plane is %0.2f to %0.2f\n',...
	lums(1), lums(end));

% Now we have the coordinates we desire, make a picture of the
% isoluminant plane.
[X,Y] = meshgrid(0:imageSize-1,0:imageSize-1);
X = X-(imageSize-1)/2; Y = Y-(imageSize-1)/2;
X = X/max(abs(X(:))); Y = Y/max(abs(Y(:)));
XVec = reshape(X,1,imageSize^2); YVec = reshape(Y,1,imageSize^2);
imageLMS = bgLMS*ones(size(XVec))+rgConeInc*XVec+sConeInc*YVec;
[imageRGB,badIndex] = SensorToSettings(calLMS,imageLMS);
bgRGB = SensorToSettings(calLMS,bgLMS);
imageRGB(:,find(badIndex == 1)) = bgRGB(:,ones(size(find(badIndex == 1))));
rPlane = reshape(imageRGB(1,:),imageSize,imageSize);
gPlane = reshape(imageRGB(2,:),imageSize,imageSize);
bPlane = reshape(imageRGB(3,:),imageSize,imageSize);
theImage = zeros(imageSize,imageSize,3);
theImage(:,:,1) = rPlane;
theImage(:,:,2) = gPlane;
theImage(:,:,3) = bPlane;

% Show the image for illustrative purposes
figure; clf; image(theImage);
