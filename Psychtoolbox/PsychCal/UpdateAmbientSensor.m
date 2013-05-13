function cal = UpdateAmbientSensor(cal,new_ambient_sensor,ADD)
% cal = UpdateAmbientSensor(cal,new_ambient_sensor,[ADD])
%
% Update the ambient light used in the conversions.  The
% value for new_ambient_sensor should be passed in the
% same units as defined by T_sensor in the call to
% SetColorSpace.
%
%ÊIf flag ADD is true, passed ambient is added to current
% value.  Otherwise passed value replaces current value.
% ADD is false if not passed.  Use caution when setting ADD
% true -- if the ambient is changing during the experiment
% you typically don't want to keep adding multiple times.
%
% If instead you want to update in the measurement units,
% call UpdateAmbient instead.
%
% 7/7/98    dhb		   Wrote it.
% 4/5/02    dhb, ly  Update for new interface.  Internal names not changed.
% 5/2/02    dhb, kr  Add ADD flag.

% Primitive dimension check
if (size(new_ambient_sensor,1) ~= size(cal.ambient_linear,1) || ...
  size(new_ambient_sensor,2) ~= size(cal.ambient_linear,2) )
	error('Old and new ambient specifications are not of same dimension');
end

% Set default on optional argument.
if (nargin < 3 || isempty(ADD))
	ADD = 0;
end

% Update the structure
if (~ADD)
	cal.ambient_linear = new_ambient_sensor;
else
	cal.ambient_linear = cal.ambient_linear + new_ambient_sensor;
end


