function cal = SetGammaMethod(cal,gammaMode,precision)
% cal = SetGammaMethod(cal,gammaMode,[precision])
%
% Set up the gamma correction mode to be used.  Options
% are:
%   gammaMode == 0 - search table using linear interpolation via interp1.
%   gammaMode == 1 - inverse table lookup.  Fast but less accurate.
%   gammaMode == 2 - exhaustive search
%
% If gammaMode == 1, then you may specify the precision of the
% inverse table.  The default is 1000 levels.
%
% See also RefitCalGamma, CalibrateFitGamma, GamutToSettings
%
% 8/4/96  dhb  Wrote it.
% 8/21/97 dhb  Update for structure.
% 3/12/98 dhb  Change name to SetGammaCorrectMode
% 5/26/12 dhb  Add real exhaustive search mode (2). 

% Check that the needed data is available. 
gammaTable = cal.gammaTable;
gammaInput = cal.gammaInput;
if isempty(gammaTable)
	error('Calibration structure does not contain gamma data');
end

% Do the right thing depending on mode.
if gammaMode == 0
	cal.gammaMode = gammaMode;
	return;
elseif gammaMode == 1
	if nargin == 2
		precision = 1000;
	end
	iGammaTable = InvertGammaTable(gammaInput,gammaTable,precision);
	cal.gammaMode = gammaMode;
	cal.iGammaTable = iGammaTable;
elseif gammaMode == 2
    cal.gammaMode = gammaMode;
else
  error('Requested gamma inversion mode %g is not yet implemented', gammaMode);
end
	
