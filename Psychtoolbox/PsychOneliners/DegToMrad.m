function angleInMradians = DegToMrad(angleInDegrees)
% angleInMradians = DegToMrad(angleInDegrees)
%
% Convert between degrees and milliradians(mrad).
%
% 2/20/13  dhb  Wrote it.
% 9/26/13  dhb  Remove dependene on degtorad, as this is in a toolbox 
%               that not everyone has.

angleInRadians = (pi/180)*angleInDegrees;
angleInMradians = 1000*angleInRadians;
