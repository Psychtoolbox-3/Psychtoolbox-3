function lab = XYZToLab(xyz,whiteXYZ)
% lab = XYZToLab(xyz,whiteXYZ)
%
% xyz is a 3 x N matrix with xyz in columns
% whiteXYZ is a length 3 vector with the white point.
% lab is a 3 x N matrix of [Lstar;astar;bstar]
%
% Formulae are taken from Wyszecki and Stiles, page 167.
%
% xx/xx/xx    baw  Created.
% xx/xx/xx    dhb  Made compatible with version 3.5
% 10/10/93    dhb  Changed name to XYZToLab.
% 5/9/02      dhb  Improved help

% Check xyz dimensions
[m,n] = size(xyz);
if ( m ~= 3 )
  error('Array xyz must have three rows')
end

% Check white point dimensions
[m,n] = size(whiteXYZ);
if ( m ~= 3 || n ~= 1)
  error('Array white is not a three vector')
end

% Separate out the compontents
X = xyz(1,:); Y = xyz(2,:); Z = xyz(3,:);
Xn = whiteXYZ(1); Yn = whiteXYZ(2); Zn = whiteXYZ(3);

% Allocate space
[m,n] = size(xyz);
lab = zeros(3,n);
fX = zeros(1,n);
fY = zeros(1,n);
fZ = zeros(1,n);

% Calculate fX
lX = find( (X/Xn) < 0.008856 );
bX = find( (X/Xn) >= 0.008856);
if (length(bX) ~= 0) 
 fX(bX) = (X(bX)/Xn) .^(1/3);
end
if (length(lX) ~= 0)
  fX(lX) = 7.787*(X(lX)/Xn) + 16/116;
end

% Calculate L
lY = find( (Y/Yn) < 0.008856 );
bY = find( (Y/Yn) >= 0.008856);
if ( length(bY) ~= 0 )
 lab(1,bY) = 116*(Y(bY)/Yn).^(1/3) - 16;
 fY(bY) = (Y(bY)/Yn) .^(1/3);
end
if ( length(lY) ~= 0 )
 lab(1,lY) = 903.3 * (Y(lY)/Yn);
 fY(lY) = 7.787*(Y(lY)/Yn) + 16/116;
end  

% Calculate fZ
lZ = find( (Z/Zn) < 0.008856 );
bZ = find( (Z/Zn) >= 0.008856);
if (length(bZ) ~= 0) 
 fZ(bZ) = (Z(bZ)/Zn) .^(1/3);
end
if (length(lZ) ~= 0)
  fZ(lZ) = 7.787*(Z(lZ)/Zn) + 16/116;
end

% Compute a and b
lab(2,:) = 500.0*(fX - fY);
lab(3,:) = 200.0*(fY - fZ);
