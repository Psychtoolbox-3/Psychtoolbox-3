% LabLuvTest
% 
% Test routines that convert to CIELAB and CIELUV.
%
% 7/11/03  dhb  Added comments, change name.

% Set some XYZ values
XYZ = [100 100 100 ; 100 75 125]';
XYZ = [XYZ , XYZ/10000000];

% Set a white point
white = [200 200 200]';

% Convert to Luv
Luv = XYZToLuv(XYZ,white);
XYZ1 = LuvToXYZ(Luv,white);

% Convert to Lab
Lab = XYZToLab(XYZ,white);
XYZ2 = LabToXYZ(Lab,white);

% Check that inversion is correct
% up to roundoff error
fix( (XYZ-XYZ1)/1e-12 )
fix( (XYZ-XYZ2)/1e-12 )
