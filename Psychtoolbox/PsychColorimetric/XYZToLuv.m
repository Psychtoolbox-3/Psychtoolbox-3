function luv = XYZToLuv(xyz,whiteXYZ)
% luv = XYZToLuv(xyz,whiteXYZ)
%
% xyz is a 3 x N matrix with xyz in columns
% whiteXYZ is a 3 vector of the white point
% luv is a 3 x N matrix with L*u*v* in the columns
%
% Formulae are taken from Wyszecki and Stiles, page 167.
%
% xx/xx/xx    baw  Created.
% xx/xx/xx    dhb  Made compatible with version 3.5
% 10/10/93    dhb  Changed name to XYZToLuv
% 5/9/02      dhb  Improved help.

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

% Compute u and v for white point
uw = (4.0 * Xn) / (Xn + 15.0*Yn + 3.0*Zn);
vw = (9.0 * Yn) / (Xn + 15.0*Yn + 3.0*Zn);

% Allocate space
[m,n] = size(xyz);
luv = zeros(m,n);

% Compute L
lY = find( (Y/Yn) < 0.008856 );
bY = find( (Y/Yn) >= 0.008856);
if ( length(bY) ~= 0 )
 luv(1,bY) = 116*(Y(bY)/Yn).^(1/3) - 16;
end
if ( length(lY) ~= 0 )
 luv(1,lY) = 903.3 * (Y(lY)/Yn);
end  

% Compute u and v
uv = XYZTouv(xyz);

% Compute u* and v*
luv(2,:) = 13 * luv(1,:) .* (uv(1,:) - uw);
luv(3,:) = 13 * luv(1,:) .* (uv(2,:) - vw);
