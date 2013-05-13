function out = Rot3d(in,n,dim)
% out = Rot3d(in,n,dim)
% rotates input matrix around x-, y- of z-axis
% IN is a vector, a 2D or a 3D matrix (all three can be rotated three
% dimensionally)
%
% rotation can only be performed by multiples of 90 degrees, N is the
% number of 90 degree steps by which the matrix will be rotated (clockwise,
% same as rot90())
% DIM indicates the axis to rotate about:
%   1=y; 2=x; 3=z
%
% Rot3d([1 2 3 4; 3 5 6 7],1,1)
%   ans(:,:,1) =
%        4
%        7
%   ans(:,:,2) =
%        3
%        6
%   ans(:,:,3) =
%        2
%        5
%   ans(:,:,4) =
%        1
%        3

% DN 2008

psychassert(nargin==3,'function requires three inputs');
psychassert(ndims(in)<=3,'input cannot have more than three dimensions');
psychassert(dim>0 && dim<=3,'invalid rotation axis specified (%d)\nuse 1=y; 2=x; 3=z',dim);

ss          = AltSize(in,[1 2 3]);
[y,x,z]     = meshgrid(1:ss(1),1:ss(2),1:ss(3));
yi          = y(:);
xi          = x(:);
zi          = z(:);

incoords    = [xi.'; yi.'; zi.']; % change to a 3xN matrix with [x; y; z]
rot         = n*90;

switch dim
    case 1
        oc = Roty(incoords,rot);
    case 2
        oc = Rotx(incoords,rot);
    case 3
        oc = Rotz(incoords,rot);
end

% make all coordinates positive
rijplus     = (max(abs(oc),[],2)+1) .* max(oc<=0,[],2);

yo          = oc(2,:)+rijplus(2);
xo          = oc(1,:)+rijplus(1);
zo          = oc(3,:)+rijplus(3);

% use original and rotated indices to create output
out         = zeros(max(yo),max(xo),max(zo));
ii          = sub2ind(ss,yi,xi,zi);
oi          = sub2ind(size(out),yo,xo,zo);

out(oi)     = in(ii);


% subfunctions (multiplication with rotation matrices)
function [XYZ] = Rotx(XYZ,a)

Rx  = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
XYZ = Rx * XYZ;

function [XYZ] = Roty(XYZ,b)

Ry  = [cosd(b) 0 sind(b); 0 1 0; -sind(b) 0 cosd(b)];
XYZ = Ry * XYZ;

function [XYZ] = Rotz(XYZ,g)

Rz  = [cosd(g) -sind(g) 0; sind(g) cosd(g) 0; 0 0 1];
XYZ = Rz * XYZ;
