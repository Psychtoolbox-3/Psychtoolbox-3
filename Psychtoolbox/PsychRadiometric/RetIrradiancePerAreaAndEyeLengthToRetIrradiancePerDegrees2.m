function retIrradiance_PerDegrees2 = RetIrradiancePerAreaAndEyeLengthToRetIrradiancePerDegrees2(retIrradiance_PerArea,eyeLength)
% retIrradiance_PerDegrees2 = RetIrradiancePerAreaAndEyeLengthToRetIrradiancePerDegrees2(retIrradiance_PerArea,eyeLength)
%
% Convert retinal irradiance measured in units of Y/x^2 to units of
% Y/deg^2, where x is a unit of distance (m, cm, mm, um, etc.) and
% Y is a measure of light amount (Watts, Joules, quanta/sec, quanta, etc.)
%
% Eye length should be passed in units of x.
%
% The conversion assumes that we are in the small angle regime, where
% degrees are essentially linear with retinal extent.
%
% See also: PsychRadiometric, RetIrradiancePerDegrees2AndEyeLengthToRetIrradiancePerArea.
%
% 6/23/13  dhb  Wrote it.

% Convert x to degrees.  The routine RetinalMMToDegrees does not
% actually care whether the input is in mm, it just needs its
% two arguments to be in the same units.
%
% We use DegreesToRetinalMM and invert it, rather than RetinalMMToDegrees.
% This is because with the former, we can force the number of degrees to be small
% and get a factor valid in the small angle range.  If we do it the
% other way, when the units of distance are large relative to the 
% eye length, weird things can happen.
xPerDegree = DegreesToRetinalMM(1,eyeLength);
degreesPerX = 1/xPerDegree;
degrees2PerArea = degreesPerX^2;

retIrradiance_PerDegrees2 = retIrradiance_PerArea/degrees2PerArea;

