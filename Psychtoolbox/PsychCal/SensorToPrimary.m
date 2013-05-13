function [primary] = SensorToPrimary(cal,sensor)
% [primary] = SensorToPrimary(cal,sensor)
%
% Convert from sensor color space coordinates to primary
% coordinates.
%
% This depends on the standard calibration globals.
%
% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow device characterization dimensions to exceed
%                  linear settings dimensions.
% 11/11/93   dhb   Update for new calData routines.
% 11/17/93   dhb   Newer calData routines.
% 8/4/96     dhb   Update for stuff bag routines.
% 8/21/97    dhb   Update for structures.
% 4/5/02     dhb, ly  New calling convention.  Internal naming not changed fully.


% Get size
[m,n] = size(sensor);

% Get necessary calibration data
M_linear_device = cal.M_linear_device;
ambient_linear = cal.ambient_linear;
if (isempty(M_linear_device) || isempty(ambient_linear))
	error('SetSensorColorSpace has not been called on calibration structure');
end

% Ambient corrections
[ma,na] = size(ambient_linear);
if (m ~= ma || na ~= 1)
  error('Incorrect dimensions for ambient');
end
sensora = sensor-ambient_linear*ones(1,n);

% Color space conversion
[mm,nm] = size(M_linear_device);
if (m > nm)
  error ('Incorrect dimensions for M_linear_device');
end
primary = M_linear_device*sensora;
