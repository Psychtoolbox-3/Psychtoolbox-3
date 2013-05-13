function angle = MunsellHueToAngle(H1,H2)
% angle = MunsellHueToAngle(H1,H2)
%
% Convert a standard Munsell hue designation to an angle
%
% The hue designation consists of a leading number
% (e.g. 2.5, 5, 7.5, 10) and a symbolic hue name.
%
% We convert these to angle using the convention specified
% on page 509 of Wyszecki and Stiles, 2cd edition.  Angle
% is in degrees increasing counterclockwise form the positive
% x-axis in the diagram.
%
% See also: MunsellAngleToHue
%
% dhb, ijk  Wrote it.

% Get base angle information from sybmolic hue name
switch (H2)
    case 'Y'
        index = 0;
    case 'YR'
        index = 1;
    case 'R'
        index = 2;
    case 'RP'
        index = 3;
    case 'P'
        index = 4;
    case 'PB'
        index = 5;
    case 'B'
        index = 6;
    case 'BG'
        index = 7;
    case 'G'
        index = 8;
    case 'GY'
        index = 9;
    otherwise
        error(sprintf('Illegal Munsell hue string passed: %s',H2));
end
angle10 = index*360/10;
angle0 = (index+1)*360/10;

% Now we need to interpolate between angle10 and angle0, depending on the
% number passed in H1.
if (H1 < 0 || H1 > 10)
    error(sprintf('Illegal Munsell hue number passed: %s',H1));
end
angle = interp1([10 0],[angle10 angle0],H1,'linear');
