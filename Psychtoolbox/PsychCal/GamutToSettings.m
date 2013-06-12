function settings = GamutToSettings(cal, gamut)
% settings = GamutToSettings(cal, gamut)
%
% Find the best device settings to produce
% the passed linear device coordinates.
% 
% The passed coordinates should be in the range [0,1].
% The returned settings also run from [0,1], but after
% inversion of the device's gamma measurements.
%
% The returned argument values is what you actually should
% get after quantization error.
%
% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow gamma table dimensions to exceed device settings.
% 11/11/93   dhb   Update for new calData routines.
% 8/4/96     dhb   Update for stuff bag routines.
% 8/21/97    dhb   Update for structures.
% 4/13/02	 awi   Replaced SettingsToDevice with new name SettingsToPrimary.
% 11/16/06   dhb   Adjust for [0,1] world.  Involves changing what's passed
%                  in and out
% 5/26/12    dhb   Add gammaMode == 2 case.

% Error checking
if isempty(cal.gammaTable)
	error('No gamma table present in calibration structure');
end
if isempty(cal.gammaMode)
	error('SetGammaMethod has not been called on calibration structure');
end

if cal.gammaMode == 0
	settings = GamutToSettingsSch(cal.gammaInput, cal.gammaTable, gamut);
elseif cal.gammaMode == 1
	%iGammaTable = cal.iGammaTable;
	if isempty(cal.iGammaTable)
		error('Inverse gamma table not present for gammaMode == 1');
	end
	settings = GamutToSettingsTbl(cal.iGammaTable, gamut);
elseif cal.gammaMode == 2
    settings = GamutToSettingsExhaust(cal.gammaInput, cal.gammaTable, gamut);
else
	error('Requested gamma inversion mode %g is not yet implemented', cal.gammaMode);
end
