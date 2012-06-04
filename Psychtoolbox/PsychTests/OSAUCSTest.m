% OSAUCSTest
%
% Test the OSA UCS conversion routines.
%
% See XYZToLjg, LjgToXYZ.
%
% 3/27/01  dhb  Wrote it.

% Establish some pairs.  These from Wyszecki and Stiles, 2cd edition,
% Table I(6.6.4).  Should match the same table in MacAdam (1978).
xyY = [[0.2588 .2169 3.22]'  [0.5132 0.3159 9.49]' [0.4420 0.4850 54.81]'];
Ljg = [[-7 -3 -1]' [-4 2 -8]' [3 11 1]'];
XYZ = xyYToXYZ(xyY);

xyY1 = [[0.4982 .4237 29.29]' [0.4691 0.4427 30.78]' [0.2618 0.2405 22.61]'];
Ljg1 = [[0 8 -4]' [0 8 -2]' [-1 -5 -1]'];
XYZ1 = xyYToXYZ(xyY1);

xyY2 = [0.3138 0.33098 20.1]';
Ljg2 = [-2 0 0]';
XYZ2 = xyYToXYZ(xyY2);

% Forward conversion.  Can check that obtained values
% are close.
Ljg_got = XYZToLjg(XYZ);
Ljg_got1 = XYZToLjg(XYZ1);
Ljg_got2 = XYZToLjg(XYZ2);

% Backward conversion.  Can check that obtained values
% are close, and also that Ljg are as desired.
XYZ_got = LjgToXYZ(Ljg);
XYZToLjg(XYZ_got); 
XYZ_got1 = LjgToXYZ(Ljg1);
XYZToLjg(XYZ_got1);
XYZ_got2 = LjgToXYZ(Ljg2);
XYZToLjg(XYZ_got2);

% Dump a difference
Ljg_got-Ljg
