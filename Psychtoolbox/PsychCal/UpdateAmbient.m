function cal = UpdateAmbient(cal,newP_ambient,ADD)
% cal = UpdateAmbient(cal,newP_ambient,[ADD])
%
% Update the ambient light used in the conversions.  The
% value of P_ambient in the structure is replaced with the
% passed value and the computed quantities that depend on
% it are updated.
%
% The ambient must be specified in the same measurement units as it
% was in the cal at the initial call to SetColorSpace.  If different
% units are desired, all ambient fields in the structure must be updated
% and SetColorSpace called again.  I have never wanted to do this,
% so I haven't written a separate routine.
%
% It is sometimes useful, however, to update the ambient in the
% linear color space defined by the call to SetColorSpace.  To
% do this, use UpdateAmbientLinear.
%
%ÊIf flag ADD is true, passed ambient is added to current
% value.  Otherwise passed value replaces current value.
% ADD is false if not passed.  Use caution when setting ADD
% true -- if the ambient is changing during the experiment
% you typically don't want to keep adding multiple times.
%
% 11/17/93  dhb		Wrote it.
% 8/4/96    dhb   Updated for modern scheme.
% 8/21/97   dhb   Update for structures.
% 3/2/98		dhb		Fix bug in checks introduce 8/21/97, pointed out by dgp.
% 3/10/98		dhb		Change T_ to P_.
% 10/26/99  dhb, mdr  Fix bug in checks. There was also a variable name
%									glitch.  I don't think this could have worked the way
%									it was.  Perhaps no one calls it.
% 5/2/02    dhb, kr  Add ADD flag.

% Set default on optional argument.
if (nargin < 3 || isempty(ADD))
	ADD = 0;
end

% Check that passed data are compatible
oldP_ambient = cal.P_ambient;
T_ambient = cal.T_ambient;
S_ambient = cal.S_ambient;
if (isempty(oldP_ambient) || isempty(T_ambient) || isempty(S_ambient))
	error('Calibration structure does not contain ambient data');
end
[nOld,mOld] = size(oldP_ambient);
[nNew,mNew] = size(newP_ambient);
if (nOld ~= nNew || mOld ~= mNew)
	error('Old and new ambient specifications are not in same units');
end

% Update
if (~ADD)
	cal.P_ambient = newP_ambient;
else
	cal.P_ambient = cal.P_ambient + newP_ambient;
end

% Get conversion matrix and convert
M_ambient_linear = cal.M_ambient_linear;
if (isempty(M_ambient_linear))
	error('SetColorSpace has not been called on this calibration structure');
end
ambient_linear = M_ambient_linear*cal.P_ambient;
cal.ambient_linear = ambient_linear;


