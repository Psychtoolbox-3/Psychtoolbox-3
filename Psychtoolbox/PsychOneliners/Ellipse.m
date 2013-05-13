function bool = Ellipse(a,b,horpow,verpow)
% Ellipse(a) creates a Circle with
% diameter == ceil(2*a)
%
% Ellipse(a,b) creates an Ellipse with
% horizontal axis == ceil(2*a) and vertical axis == ceil(2*b)
%   
% Ellipse(a,b,power) generates a superEllipse according to the
% geometric formula (x./a).^power + (y./b).^power < 1
%
% Ellipse(a,b,horpow,verpow) generates a generalized superEllipse according
% to the geometric formula(x./a).^horpow + (y./b).^verpow < 1
%
% For more info on superEllipses, see
%   http://en.wikipedia.org/wiki/SuperEllipse
%
% Ellipse returns a (tighly-fitting) boolean matrix which is true for all
% points on the surface of the Ellipse and false elsewhere

% DN 2008
% DN 2009-02-02 Updated to do Circles and input argument handling more
%               efficiently
% DN 2011-08-31 Output wasn't always of right size (ceil(2*input))

error(nargchk(1, 4, nargin, 'struct'));

if nargin < 2
    b = a;
end
if nargin < 3
    horpow = 2;
end
if nargin < 4
    verpow = horpow;
end

[x,y]   = meshgrid(linspace(-a,a,ceil(2*a)+2),linspace(-b,b,ceil(2*b)+2));

bool    = abs(x./a).^horpow + abs(y./b).^verpow  < 1;

% return in a tight-fitting matrix
cropcoords  = CropBlackEdges(bool);
bool        = bool(cropcoords(3):cropcoords(4),cropcoords(1):cropcoords(2));
