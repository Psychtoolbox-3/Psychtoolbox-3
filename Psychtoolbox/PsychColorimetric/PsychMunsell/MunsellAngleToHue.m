function [H1,H2] = MunsellAngleToHue(angle)
% [H1,H2] = MunsellAngleToHue(angle)
%
% Invert MunsellHueToAngle.
%
% See also: MunsellHueToAngle
%
% dhb, ijk  Wrote it.

% Put angle into range 0 <= angle < 360
while (angle >= 360)
    angle = angle-360;
end
while (angle < 0)
    angle = angle+360;
end

% Figure out which of the 10 symbolic hue names we've got
if (0 <= angle && angle < 360/10)
    H2 = 'Y';
    remainder = angle - 0;
elseif (360/10 <= angle && angle < 2*360/10)
    H2 = 'YR';
    remainder = angle - 360/10;
elseif (2*360/10 <= angle && angle < 3*360/10)
    H2 = 'R';
    remainder = angle - 2*360/10;
elseif (3*360/10 <= angle && angle < 4*360/10)
    H2 = 'RP';
    remainder = angle - 3*360/10;
elseif (4*360/10 <= angle && angle < 5*360/10)
    H2 = 'P';
    remainder = angle - 4*360/10;
elseif (5*360/10 <= angle && angle < 6*360/10)
    H2 = 'PB';
    remainder = angle - 5*360/10;
elseif (6*360/10 <= angle && angle < 7*360/10)
    H2 = 'B';
    remainder = angle - 6*360/10;
elseif (7*360/10 <= angle && angle < 8*360/10)
    H2 = 'BG';
    remainder = angle - 7*360/10;
elseif (8*360/10 <= angle && angle < 9*360/10)
    H2 = 'G';
    remainder = angle - 8*360/10;
elseif (9*360/10 <= angle && angle < 10*360/10)
    H2 = 'GY';
    remainder = angle - 9*360/10;
else
    error('Logic error in routine, should not be possible to be here');
end
    
% Now we need to interpolate the remainder to get H1
% number passed in H1.
if (remainder < 0 || remainder >= 360/10)
    error(sprintf('Illegal Munsell hue number passed: %s',H1));
end
H1 = interp1([0 360/10],[10 0],remainder,'linear');
