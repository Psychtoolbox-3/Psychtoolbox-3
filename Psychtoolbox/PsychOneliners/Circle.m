function bool = Circle(rad)
% cirbool = Circle(rad)
%
% Circle creates a Circle with diameter == ceil(2*rad)
% Ellipse returns a (tighly-fitting) boolean matrix which is true for all
% points on the surface of the Ellipse and false elsewhere

% DN 2008
% DN 2009-02-02 Turned into proxy for Ellipse

bool = Ellipse(rad);
