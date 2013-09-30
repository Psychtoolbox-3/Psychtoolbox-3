function angleInDegrees = MradToDeg(angleInMradians)
% angleInDegrees = MradToDeg(angleInMradians)
%
% Convert between milliradians(mrad) and degrees.
%
% 2/20/13  dhb  Wrote it.
% 9/26/13  dhb  Remove dependene on radtodeg, as this is in a toolbox 
%               that not everyone has.

angleInRadians = angleInMradians/1000;
angleInDegrees = (180/pi)*angleInRadians;
