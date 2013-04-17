function angleInMradians = DegToMrad(angleInDegrees)
% angleInMradians = DegToMrad(angleInDegrees)
%
% Convert between degrees and milliradians(mrad).
%
% 2/20/13  dhb  Wrote it.

angleInMradians = 1000*degtorad(angleInDegrees);
