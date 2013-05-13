function cal = DropCalBits(cal,whichScreen,forceBits)
% cal = DropCalBits(cal,whichScreen,[forceBit])
%
% Drops the bitdepth of a calibration file if
% necessary.  Useful for running programs
% transparently on 8 and 10 bit hardware.
%
% If arg forceBits is passed, it is used as
% the current hardware depth.  Otherwise the
% reported DACBits of whichScreen is used.
%
% This code assumes calibration was done at
% equally spaced levels in RGB settings, as is
% the case with our calibration routines.  May
% not generalize, and I haven't worried about
% the roundoff errors.  Certainly OK for basic
% use.
%
% 2/13/05		dhb		Wrote it.

% Get hardware dac level.  Note that the application
% code should use LoadClut, not SetClut, to access
% full bit depth.
if (nargin > 2 && ~isempty(forceBits))
	hardwareBits = forceBits;
else
	hardwareBits = Screen(whichScreen,'Preference','DACBits');
end

% Force calibration down to 8 bits, which is how we plan to use it.
% Simply refit raw data at correct number of input levels.  
if (cal.describe.dacsize > hardwareBits)
	cal.describe.dacsize = hardwareBits;
	nInputLevels = 2^cal.describe.dacsize;
	cal.rawdata.rawGammaInput = round(linspace(nInputLevels/cal.describe.nMeas,nInputLevels-1,cal.describe.nMeas))';
	cal = CalibrateFitGamma(cal);
elseif (cal.describe.dacsize < hardwareBits)
	error('Current hardware has greater bit depth than at calibration.');
end
	
