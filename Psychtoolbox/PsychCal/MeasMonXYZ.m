function XYZ = MeasMonXYZ(window,settings,syncMode,whichMeterType)
% XYZ = MeasMonXYZ(window,settings,[syncMode],[whichMeterType])
%
% Measure the XYZ of a series of monitor settings.
%
% This routine is specific to go with CalibrateMon,
% as it depends on the action of SetMon. 
%
% If whichMeterType is passed and set to 0, then the routine
% returns random spectra.  This is useful for testing when
% you don't have a meter.
%
% Other valid types:
%  1 - Use PR650 (default)
%
% 10/26/93  dhb  	Wrote it based on ccc code.
% 11/12/93  dhb  	Modified to use SetColor.
%	6/23/94		ccc		Modified it from MeasMonSpd.m for
%									the purpose of measuring XYZ
%	8/9/94		dhb		Added code to go into sync mode
%									And then commented it out.
% 8/11/94		dhb		Sync mode back in
% 8/15/94		dhb		Sync mode as argument.
% 4/12/97   dhb   New toolbox compatibility, take window and bits args.
% 8/26/97		dhb, pbe Add noMeterAvail option.
% 4/7/99    dhb   Add argument for radius board.  Compact default arg code.
% 8/14/00   dhb   Call to CMETER('SetParams') conditional on OS9.
% 8/20/00   dhb   Remove bits arg from call to SetColor.
% 8/21/00   dhb   Remove dependence on RADIUS flag.  This is now handled inside of SetColor.
%	          dhb   Change calling conventions to remove unused args.
% 9/14/00   dhb   Sync mode no longer used.  Arg passed for backwards compatibility.
% 2/27/02   dhb, ly  Pass whichMeterType rather than noMeterAvail.

% Check args and make sure window is passed right.
usageStr = 'XYZ = MeasMonXYZ(window,settings,[syncMode],[whichMeterType])';
if (nargin < 2 || nargin > 4 || nargout > 1)
	error(usageStr);
end
if (size(window,1) ~= 1 || size(window,2) ~= 1)
	error(usageStr);
end

% Set defaults
defaultSync = 0;
defaultWhichMeterType = 1;

% Check args and set defaults
if (nargin < 4 || isempty(whichMeterType))
	whichMeterType = defaultWhichMeterType;
end
if (nargin < 3 || isempty(syncMode))
	syncMode = defaultSync;
end

[null,nMeas] = size(settings);
XYZ = zeros(3,nMeas);
for i=1:nMeas
	% Set color
  SetColor(window,1,settings(:,i));
	
	% Make the measurement
	switch (whichMeterType)
		case 0,
			XYZ(:,i) = sum(settings(:,i)*ones(3,1);
			WaitSecs(0.1);
		case 1,
		  XYZ(:,i) = MeasXYZ;
		otherwise,
			error('Invalid meter type specified');
	end
end

