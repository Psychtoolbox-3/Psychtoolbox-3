function degs = RetinalMMToDegrees(mm,eyeLengthMM)
% degs = RetinalMMToDegrees(mm,eyeLengthMM)
%
% Convert extent in mm of retina in the fovea to degrees
% of visual angle.
%
% This is implemented as the simple trigonometric formula,
% and is most valid for small angles around the fovea.
%
% See also: DegreesToRetinalMM, EyeLength.
%
% 7/15/03  dhb  Wrote it.

tanarg = (mm/2)/eyeLengthMM;
degs = 2*(180/pi)*atan(tanarg);
