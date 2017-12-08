% Psychtoolbox:PsychColorimetric.
%
% Colorimetric and spectral calculations.  See also Psychtoolbox:PsychColorimetricData.
% And also see the closely related Psychtoolbox:PsychRadiometric.  Sometimes it is
% not entirely clear whether a routine should be classified as colorimetric or
% radiometric. Apologies if our intuitions don't match yours
%
% help Psychtoolbox              % For an overview, triple-click me & hit enter.
% help PsychColorimetricData     % For colorimetric calculations, triple-click me & hit enter.
% help PsychColorimetricMatFiles % For list of data .mat files, triple-click me & hit enter.
% help PsychRadiometric          % For list of data .mat files, triple-click me & hit enter.
%
%   AbsorbanceToAbsorbtance - Obsolete.  Use AbsorbanceToAbsorptance (spelling fixed).
%   AbsorbanceToAbsorptance - Convert absorbance to absorptance spectrum.
%   AbsorbtanceToAbsorbance - Obsolete.  Use AbsorptanceToAbsorbance (spelling fixed).
%   AbsorptanceToAbsorbance - Convert absorptance to absorbance spectrum.
%   CheckWls            - Check consistency of two wavelength descriptions.
%   ComputeDE           - Compute vector length between matrix columns.
%   ComputeDE2000_Lab   - Compute deltaE values (CIEDE200) for color pairs given in L*a*b* coordinates.
%   ComputeDKL_M        - Compute transformation matrix for DKL color space.
%   ComputeAxialDensity - Compute axial optical density from specific density and path length.
%   ComputePhotopigmentBleaching - Compute fraction of photopigment bleached.
%   ConeIncToDKL        - Convert from cone increments to DKL.
%   ContrastToExcitation - Convert contrast to excitation coordinate.
%   ContrastToIncrement - Convert contrast to incremental coordinates.
%   DKLToConeInc        - Convert from DKL to cone increments.
%   EffectiveTrolandsFromLum - Compute effective trolands from luminance.
%   EnergyToCones       - Convert monochromatic energy to cone excitations.
%   ExcitationToContrast - Convert excitation coordinates to contrast.
%   FindCmfPeaks        - Find the peak wavelength and values of a set of color matching functions.
%   FindCovLinMod       - Find linear model from covariance matrix.
%   FindLinMod          - Find linear model for an ensemble of vectors.
%   FindModelWeights    - Find the weights with respect to a linear model.
%   FToRatio            - Subroutine for Lab/Luv conversions.
%   GenerateBlackBody   - Generate spectra for black body radiators.
%   GenerateCIEDay      - Generate CIE daylights.
%   IncrementToContrast - Convert incremental coordinates to contrast.
%   IsomerizationsFromAbsorbptions - Compute isomerization rate from absorption rate.
%   LabToXYZ            - Convert from Lab to XYZ.
%   LjgToXYZ            - Convert from OSA UCS Ljg to XYZ (10 degree).
%   LjgToXYZFun         - Error function for LjgToXYZ numerical optimization.
%   LumToRadiance       - Get spectral radiance from luminance and relative spectrum of source.
%   LumToTrolands       - Convert luminance (cd/m2) to trolands.
%   LuvToXYZ            - Convert from Luv to XYZ.
%   LMSToMacBoyn        - Convert from cones to MacLeod-Boynton chromaticity.
%   MakeFourierBasis    - Make a set of Fourier component basis functions.
%   MakeGaussBasis      - Make a set of Gaussian basis functions.
%   MakeItS             - Force wavelength sampling spec. to S format.
%   MakeItStruct        - Force wavelength sampling spec. to struct foramt.
%   MakeItWls           - Force wavelength sampling spec. to wls format.
%   MakeMonoPrimary     - Make the spd of a monochromatic primary.
%   MakeOrtho           - Make a set of basis functions orthonormal.
%   MakeUnitLength      - Make the columns of a matrix have unit length.
%   MonoImageToSRGB     - Convert a monochrome image to an sRGB color image at passed chromaticity.
%   M_PToP              - Conversion matrix from source/dest. primaries.
%   M_PToT              - Conversion matrix from source primaries and dest. cmfs.
%   M_TToP              - Conversion matrix from source cmfs and dest. primaries.
%   M_TToT              - Conversion matrix from source/dest. cmfs.
%   PhotonAbsorptionRate - Compute photon absoroption rate.
%   PupilDiameterFromLum - Estimate pupil diameter from luminance.
%   PsychMunsell        - Munsell renotation to xyY conversion.
%   RetIrradianceToIsoRecSec - Convert retinal irradiance (power units) to iso. per receptor per second.
%   SampleCircle        - Sample points on a circle.
%   SampleSphere        - Sample points on a sphere.
%   ShiftSpectra        - Shift a spectral function along the wavelength axis.
%   SPDToLinSPD         - Find approximation within specified linear model.
%   SPDToMetSPD         - Find a metameric spd.
%   SplineCmf           - Spline color matching functions to new wavelength sampling.
%   SplineRaw           - Subroutine for other spline functions.
%   SplineSpd           - Spline a spectral power distribution to new wavelength sampling.
%   SplineSrf           - Spline a surface reflectance function to new wavelength sampling.
%   SRGBPrimaryToXYZ    - Convert between sRGB primary coordinates and XYZ.
%   SRGBGammaCorrect    - Convert between sRGB primary coordinates and 8-bit RGB values.
%   SRGBGammaUncorrect   - Convert between sRGB 8-bit RGB values and primary coordinates.
%   SToWls              - Convert S wavelength sampling spec to wls format.
%   TestLxx             - Test routine for Lab/Luv calculations.
%   TriToMetSPD         - Compute metamer from tristimulus coordinates.
%   WattsToRetIrradiance - Get absolute retinal irradiance (power units) from rel. spectrum and watts/area.
%   uvTols              - Convert between CIE u'v' and a cone based (ls) chromaticity.
%   uvToxy              - Convert between CIE u'v' and CIE xy chromaticity.
%   uvYToXYZ            - Convert between u'v'Y and XYZ.
%   WlsToS              - Convert wls wavelength sampling spec to S format. 
%   WlsToT              - Compute identity color matching matrix from wavelength sampling.
%   xyTouv              - Convert CIE xy chromaticity to CIE u'v' chromaticity.
%   xyYToXYZ            - Convert between xyY and XYZ.
%   XYZToF              - Subroutine for Lab/Luv calculations.
%   XYZToLab            - Convert between XYZ and Lab.
%   XYZToLjg            - Convert between XYZ (10 degree) and OSA UCS Ljg.
%   XYZToLuv            - Convert between XYZ and Luv.
%   XYZToSRGBPrimary    - Convert between XYZ and sRGB primary coordinates.
%   XYZTouv             - Compute uv chromaticities from XYZ.
%   XYZTouvY            - Convert between XYZ and u'v'Y.
%   XYZToxyY            - Convert between XYZ and xyY.
  
% Copyright (c) 1996-2013 by David Brainard, Denis Pelli, & Mario Kleiner


