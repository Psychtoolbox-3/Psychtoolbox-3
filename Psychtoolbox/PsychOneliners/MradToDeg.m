function angleInDegrees = MradToDeg(angleInMradians)
% angleInDegrees = MradToDeg(angleInMradians)
%
% Convert between milliradians(mrad) and degrees.
%
% 2/20/13  dhb  Wrote it.

angleInDegrees = radtodeg(angleInMradians/1000);
