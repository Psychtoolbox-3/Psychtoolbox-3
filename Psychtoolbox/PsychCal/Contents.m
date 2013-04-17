% Psychtoolbox:PsychCal.
%
% Calibration routines for monitors, e.g., gamma correction.
% For a very simple gamma calibration of a standard monitor with a manually
% controlled photometer, use CalibrateMonitorPhotometer.
%
% For advanced calibration, with manual or automatic photometers, with
% standard monitors or high precision displays like Bits+,
% DataPixx/ViewPixx etc. use CalibrateMonSpd.
%
%
% help Psychtoolbox % For an overview, triple-click me & hit enter.
% help PsychDemos   % For demos, triple-click me & hit enter.
%
%   BasicToneMapCalFormat - Apply simple (truncation) tone mapping to XYZ data in cal format.
%   CalDataFolder       - Path to the calibration data folder "PsychCalData".
%   CalFormatToImage    - Put a stretched out image back into an image image.
%   CalibrateAmbDrvr    - Common code called by monitor calibration programs.  Script.
%   CalibrateFitGamma   - Fit the gamma function to the calibration measurements.
%   CalibrateFitLinMod  - Fit the linear model to spectral calibration data.
%   CalibrateFitYoked   - Fit a very special kind of yoked calibration data.
%   CalibrateManualDrvr - Driver for manual monitor calibration, called by CalibrateMonSpd.
%   CalibrateMonDrvr    - Common code called by monitor calibration programs.  Script.
%   CalibrateMonitorPhotometer - Simple calibration procedure for use with manually controlled photometers.
%   CalibrateMonSpd     - Run standard monitor spectral calibration.  Script.
%   CalibratePlotGamma  - Plot the gamma functions in a calibration structure.
%   CalibratePlotSpectra - Plot the device spectra in a calibration structure.
%   CalToggleBitsPlusPlus - Determine whether to calibrate with or without Bits++.
%   CompareMonCal       - Compare two calibration structures to see if they match.
%   CompareMonCalOverTime - Compare two calibrations of same monitor and see differences.  Script.
%   ContrastMatch       - Match contrast of two gratings.
%   Cross2D             - Create 2D vectors that cross combinations of passed vectors.
%   Cross3D             - Create 3D vectors that cross combinations of passed vectors.
%   CVIAdjustCalData    - Adjust of PR650 measurements with CVI measurements.  Not of general interest.
%   CylToSensor         - Convert from cylindrical to sensor coordinates.
%   DescribeMonCal      - Print a description of calibration structure.
%   DumpMonCalSpd       - Dump contents of standard spectral calibration file.  Script.
%   ExpandSettings      - Subroutine for "Acc" conversion routines.
%   FindSpectralPeaks   - Find the peaks in a spectral function.
%   FlushCalFile        - Remove old entries from calibration file.
%   GamutToSettings     - Convert from within gamut device coordinates to settings.
%   GamutToSettingsSch  - Subroutine for GamutToSettings.
%   GamutToSettingsTbl  - Subroutine for GamutToSettings.
%   GratingNull         - Match two luminances.
%   IccProfile.mex      - Read ICC profiles for monitor calibration.
%   ImageToCalFormat    - Stretch out an image to pass to calibration conversion routines.
%   InvertGammaTable    - Subroutine for GamutToSettings.
%   LoadCalFile         - Load calibration data.
%   MaximizeGamutContrast- Maximize contrast within gamut.
%   MeasMonSpd          - Set monitor and measure spd.
%   MeasMonXYZ          - Set monitor and measure XYZ.
%   MeasureDpi          - Measure monitor dpi.
%   PolarToSensor       - Convert from polar to linear coordinates.
%   PrimaryToGamut      - Force primary coordinates into gamut.
%   PrimaryToSensor     - Convert from primary to sensor coordinates.
%   PrimaryToSettings   - Convert from primary coordinates to device settings.
%   RefitCalGamma       - Refit the gamma data in a monitor calibration file.  Script.
%   RefitCalLinMod      - Refit primary spectral linear model in calibration file.  Script.
%   RefitCalYoked       - Refit linear model to match yoked [a,a,a] measurements. Script.
%   SaveCalFile         - Save calibration data.
%   SearchGammaTable     - Subroutine for GamutToSettings.
%   SensorToCyl         - Convert from sensor to cylindrical coordinates.
%   SensorToPolar       - Convert from sensor to polar coordinates.
%   SensorToPrimary     - Convert from sensor to primary coordinates.
%   SensorToSettings    - Convert from sensor coordinates to device settings.
%   SensorToSettingsAcc - Convert from sensor coordinates to settings, handle spectral shifts.
%   SetGammaMethod      - Set gamma correction method for conversions.
%   SetSensorColorSpace - Set working sensor color space for conversions.
%   SettingsToPrimary   - Convert from settings to primary coordinates.
%   SettingsToSensor    - Convert from settings to sensor coordinates.
%   SettingsToSensorAcc - Convert from settings to sensor, handle spectral shifts.
%   UpdateAmbient       - Update the ambient used for conversions, in meas. units.
%   UpdateAmbientSensor - Update the ambient used for conversions, in sensor units.

% Copyright (c) 1996-2010 by Denis Pelli & David Brainard

