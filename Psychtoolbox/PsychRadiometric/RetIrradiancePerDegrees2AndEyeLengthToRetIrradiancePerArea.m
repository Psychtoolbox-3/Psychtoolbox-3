function retIrradiance_PerArea = RetIrradiancePerDegrees2AndEyeLengthToRetIrradiancePerArea(retIrradiance_PerDegrees2,eyeLength)
% retIrradiance_PerArea = RetIrradiancePerDegrees2AndEyeLengthToRetIrradiancePerArea(retIrradiance_PerDegrees2,eyeLength)
%
% Convert retinal irradiance measured in units of Y/deg^2 to units of
% Y/x^2, where x is a unit of distance (m, cm, mm, um, etc.) and
% Y is a measure of light amount (Watts, Joules, quanta/sec, quanta, etc.)
%
% Eye length should be passed in units of x.
%
% The conversion assumes that we are in the small angle regime, where
% degrees are essentially linear with retinal extent.
%
% See also: PsychRadiometric, RetIrradiancePerAreaAndEyeLengthToRetIrradiancePerDegrees2.
%
% 6/23/13  dhb  Wrote it.

% Convert x to degrees.  The routine DegreesToRetinalMM does not
% actually care whether the input is in mm, it just needs its
% two arguments to be in the same units.
%
% We use DegreesToRetinalMM, rather than RetinalMMToDegrees followed
% by inversion.  This is because with the former, we can force the number
% of degrees to be small and get a factor valid in the small angle range.
% If we do it the other way, when the units of distance are large relative to the 
% eye length, weird things can happen.
xPerDegree = DegreesToRetinalMM(1,eyeLength);
areaPerDegrees2 = xPerDegree^2;
retIrradiance_PerArea = retIrradiance_PerDegrees2/areaPerDegrees2;

