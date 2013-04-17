% PTBAndVSETColorimetry
% 
% This script compares values calculated by PTB and vset (developed by
% Wandell and colleague) for various colorimetric functions.
%
% For information on vset, see https://github.com/wandell/vset/blob/master/README.txt.
% It contains its own implementation of many of the colorimetric computations 
% implemented in PTB.  Note that Brainard and Wandell are not completely
% independent sources.
%
% Make sure vset is on your path to run these comparisons.  The tests should
% also work if you have iset on your path instead.
%
% The checks are grouped into cells that check one thing at a time.
%
% Because DHB and BAW have too much time on their hands -
%  "A man with one watch always knows what time it is, a man with two is never sure."
%  "A man with SPM knows how to blur and average his watches so that he can't tell the time"
%
% 7/7/10  dhb  Wrote it.
% x/xx/10 baw  Checked with ISET-4.0 revision 351 (BW)
% 2/28/13 dhb  Updated to work with VSET rather than proprietary iset.
%
% TODO:  
%   There are small differences in sRGB calculations that need to be tracked down.

%% Clear and close
clear; close all;
if (exist('ISET','file'))
    fprintf('\nComparing PTB and VSET/ISET Colorimetry\n');
else
    error('Need vset (or iset) on your path to run this test');
end

%% XYZ-related colorimetry
%
% xy
fprintf('\n***** Basic XYZ *****\n');
testXYZs = [[1 2 1]' [2 1 0.5]' [1 1 1]' [0.6 2.3 4]'];
ptbxyYs = XYZToxyY(testXYZs);
ptbxys = ptbxyYs(1:2,:);
isetxys = chromaticity(testXYZs')';
if (any(abs(ptbxys-isetxys) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for XYZ to xy\n');
else
    fprintf('PTB-ISET AGREE for XYZ to xy\n');
end

%% xyY
ptbXYZs = xyYToXYZ(ptbxyYs);
if (any(abs(testXYZs-ptbXYZs) > 1e-10))
    fprintf('PTB FAILS XYZ to xyY to XYZ\n');
else
    fprintf('PTB PASSES XYZ to xyY to XYZ\n');
end
isetXYZs = xyy2xyz(ptbxyYs')';
if (any(abs(testXYZs-isetXYZs) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for xyY to XYZ\n');
else
    fprintf('PTB-ISET AGREE for xyY to XYZ\n');
end 

%% CIE uv chromaticity
ptbuvYs = XYZTouvY(testXYZs);
ptbuvs = ptbuvYs(1:2,:);
[isetus,isetvs] = XYZ2uv(testXYZs');
isetuvs = [isetus' ; isetvs'];
if (any(abs(ptbuvs-isetuvs) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for XYZ to uv\n');
    fprintf('\tI think this is because ISET implements an obsolete version of the standard\n');
else
    fprintf('PTB-ISET AGREE for XYZ to uv\n');
end

%% CIELUV
whiteXYZ = [3,4,3]';
ptbLuvs = XYZToLuv(testXYZs,whiteXYZ);
isetLuvs = xyz2luv(testXYZs',whiteXYZ')';
if (any(abs(ptbLuvs-isetLuvs) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for XYZ to Luv\n');
    fprintf('\tPresumably because the uv transformation differs.\n');
else
    fprintf('PTB-ISET AGREE for XYZ to Luv\n');
end

%% CIELAB
whiteXYZ = [3,4,3]';
ptbLabs = XYZToLab(testXYZs,whiteXYZ);
isetLabs = xyz2lab(testXYZs',whiteXYZ')';
if (any(abs(ptbLabs-isetLabs) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for XYZ to Lab\n');
else
    fprintf('PTB-ISET AGREE for XYZ to Lab\n');
end
ptbXYZCheck = LabToXYZ(ptbLabs,whiteXYZ);
isetXYZCheck = lab2xyz(isetLabs',whiteXYZ')';
if (any(abs(testXYZs-ptbXYZCheck) > 1e-10))
    fprintf('PTB FAILS XYZ to Lab to XYZ\n');
else
    fprintf('PTB PASSES XYZ to Lab to XYZ\n');
end
if (any(abs(testXYZs-isetXYZCheck) > 1e-10))
    fprintf('ISET FAILS XYZ to Lab to XYZ\n');
else
    fprintf('ISET PASSES XYZ to Lab to XYZ\n');
end

%% sRGB
%
% The iset routines seem to use a matrix that is 1/100 of the standard
% definition.  Or, at least, 1/100 of what the PTB routines use.  To
% account for this, I multiply XYZ values by 100 before passing them
% to the iset routines.
%
% The iset routines take an exponent argument for the sRGB gamma transform.
% At http://www.w3.org/Graphics/Color/sRGB this is specified as 2.4.  The
% iset routine xyz2srgb uses 2.2, with a comment that an update at
% www.srgb.com changed this from 2.4 to 2.2.  I think this is wrong,
% though.  The text I can find on the web says that the 2.4 exponent, plus
% the linear initial region, is designed to approximate a gamma of 2.2.
% That is, you want 2.4 in the formulae to approximate the 2.2 industry
% standard gamma.  Site www.srgb.com now appears gone, by the way, but all
% the other sites I find seem to be the same in this regard.
% 
% Also note that if you change the exponent in the iset sRGB formulae, you
% also should probably change the cutoff used at the low-end, where the
% sRGB standard special cases the functional form of the gamma curve.  Here
% the test is set for 2.4.
%
% Finally,the default gamma used by iset lrgb2srgb and by xzy2srgb
% currently differ, so you really want to be careful using these.  The
% inverse routine srgb2lrgb doesn't allow passing of the exponent, and it
% is hard coded as 2.4.  This causes a failure of the iset sRGB gamma
% routines to self-invert for gamma other than 2.4, and with their
% defaults.
%
% One other convention difference is that the PTB routine rounds to
% integers for the settings, while the iset routine leaves the rounding up
% to the caller.
fprintf('\n***** sRGB *****\n');

% Create some test sRGB values and convert them in the PTB framework
ptbSRGBs = [[188 188 188]' [124 218 89]' [255 149 203]' [255 3 203]'];
ptbSRGBPrimary = SRGBGammaUncorrect(ptbSRGBs);
ptbXYZs = SRGBPrimaryToXYZ(ptbSRGBPrimary);

% The ISET form takes the frame buffer values in the [0,1] regime
isetSRGBs = ptbSRGBs/255;
isetSRGBs = XW2RGBFormat(isetSRGBs',4,1);
isetXYZ   = srgb2xyz(isetSRGBs);
isetXYZs  = RGB2XWFormat(isetXYZ)';

if (any(abs(isetXYZs-ptbXYZs) > 1e-10))
    d = isetXYZs - ptbXYZs;
    fprintf('PTB-ISET DIFFERENCE for XYZ to sRGB: %f\n',max(abs(d(:))));
    d = d ./ptbXYZs;
    fprintf('PTB-ISET Percent XYZ DIFFERENCE: %f\n',max(abs(d(:))));
else
    fprintf('PTB-ISET AGREE for XYZ to sRGB\n');
end

% PTB testing of inversion
if (any(abs(XYZToSRGBPrimary(ptbXYZs)-ptbSRGBPrimary) > 1e-10))
    fprintf('PTB FAILS linear sRGB to XYZ to linear sRGB\n');
else
    fprintf('PTB PASSES linear sRGB to XYZ to linear sRGB\n');
end
if (any(abs(SRGBGammaCorrect(ptbSRGBPrimary)-ptbSRGBs) > 1e-10))
    fprintf('PTB FAILS sRGB to linear sRGB to sRGB\n');
else
    fprintf('PTB PASSES sRGB to linear sRGB to sRGB\n');
end

% Compare sRGB matrices
[nil,ptbSRGBMatrix] = XYZToSRGBPrimary([]);
isetSRGBMatrix = colorTransformMatrix('xyz2srgb')';

if (any(abs(ptbSRGBMatrix-isetSRGBMatrix) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for sRGB transform matrix\n');
else
    fprintf('PTB-ISET AGREE for sRGB transform matrix\n');
end

% XYZ -> lRGB 
% Reformat shape
ptbXYZsImage = CalFormatToImage(ptbXYZs,1,size(ptbXYZs,2));

% ISET convert 
[isetSRGBImage,isetSRGBPrimaryImage] = xyz2srgb(ptbXYZsImage);

% Reformat shape
isetSRGBs = ImageToCalFormat(isetSRGBImage); 
isetSRGBPrimary = ImageToCalFormat(isetSRGBPrimaryImage);

if (any(abs(ptbSRGBPrimary-isetSRGBPrimary) > 1e-10))
    d = ptbSRGBPrimary - isetSRGBPrimary;
    fprintf('PTB-ISET DIFFERENCE for XYZ to sRGB: %f\n',max(abs(d(:))));
    d = d ./isetSRGBPrimary;
    fprintf('PTB-ISET Percent RGB DIFFERENCE: %f\n',max(abs(d(:))));
else
    fprintf('PTB-ISET AGREE for XYZ to linear sRGB\n');
end

% ISET/PTB sRGB comparison in integer gamma corrected space
if (any(abs(round(isetSRGBs*255)-ptbSRGBs) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for XYZ to sRGB\n');
else
    fprintf('PTB-ISET AGREE for XYZ to sRGB\n');
end

% lrgb -> srgb -> lrgb in ISET
isetSRGBPrimaryCheckImage = srgb2lrgb(isetSRGBImage);
isetSRGBPrimaryCheck = ImageToCalFormat(isetSRGBPrimaryCheckImage);
if (any(abs(isetSRGBPrimaryCheck-isetSRGBPrimary) > 1e-10))
    fprintf('ISET FAILS linear sRGB to sRGB to linear sRGB\n');
else
    fprintf('ISET PASSES linear sRGB to sRGB to linear sRGB\n');
end

%% Quanta/Energy
%
% The ISET routines define c and h to more places than the PTB, so the
% agreement is only good to about 5 significant places.  Seems OK to me.
fprintf('\n***** Energy/Quanta *****\n');
load spd_D65
spdEnergyTest = spd_D65;
wlsTest = SToWls(S_D65);
testPlaces = 5;
ptbQuanta = EnergyToQuanta(wlsTest,spdEnergyTest);
isetQuanta = Energy2Quanta(wlsTest,spdEnergyTest);
if (any(abs(ptbQuanta-isetQuanta) > (10^-testPlaces)*min(ptbQuanta)))
    fprintf('PTB-ISET DIFFERENCE for energy to quanta conversion at %d significant places\n',testPlaces);
else
    fprintf('PTB-ISET AGREE for energy to quanta conversion to %d significant places\n',testPlaces);
end
if (any(abs(QuantaToEnergy(wlsTest,ptbQuanta)-spdEnergyTest) > 1e-10))
    fprintf('PTB FAILS energy to quanta to energy\n');
else
    fprintf('PTB PASSES energy to quanta to energy\n');
end
if (any(abs(Quanta2Energy(wlsTest,isetQuanta')'- spdEnergyTest) > 1e-10))
    fprintf('ISET FAILS energy to quanta to energy\n');
else
    fprintf('ISET PASSES energy to quanta to energy\n');
end

%% CIE daylights
% 
% These routines are now running in ISET and everything agrees.
fprintf('\n***** CIE Daylights *****\n');
load B_cieday
testWls = SToWls(S_cieday);
testTemp = 4987;
ptbDaySpd = GenerateCIEDay(testTemp,B_cieday);
ptbDaySpd = ptbDaySpd/max(ptbDaySpd(:));

% Iset version of normalized daylight
isetDaySpd = daylight(testWls,testTemp);
isetDaySpd = isetDaySpd/max(isetDaySpd(:));
if (any(abs(isetDaySpd-ptbDaySpd) > 1e-10))
    fprintf('PTB-ISET DIFFERENCE for daylight');
else
    fprintf('PTB-ISET AGREE for for daylight\n');
end

%% Not yet really implemented
% The code below calls some vset/iset routines that
% either don't yet have PTB counterpards, or for which
% we have not yet written the comparison code.

% Calculate correlated color temperature
isetColorTemp = cct([.31,.32]');

% Calculate spd of a daylight with some color temperature
wave = 400:700;
Tc = 6500; spd = cct2sun(wave, Tc, 'energy'); % plot(wave,spd);

% deltaE94.m
% deltaE2000.m
% deltaEab.m
% deltaEuv.m

% ieResponsivityConvert.m
% ieScotopicLuminanceFromEnergy.m
bb = blackbody(wave,6500,'watts');
chromaticity(ieXYZFromEnergy(bb',wave));
ieLuminanceFromEnergy(bb',wave);

bb = blackbody(wave,6500,'quanta');
chromaticity(ieXYZFromPhotons(bb',wave));
ieLuminanceFromPhotons(bb',wave);

