function [cal, cals, fullFilename] = LoadCalFile(filespec, whichCal, dir)
% [cal, cals, fullFilename] = LoadCalFile([filespec], [whichCal], [dir])
%
% Load calibration data from saved file in the CalData folder.
% Will search one level deep in the CalData folder if the 
% file does not exist at the top level, but skips subdirs
% called 'xOld', 'Plots', and those that begin with '.'.
%
% If no argument is given, loads from file default.mat.  If
% an integer N is passed, loads from file screenN.mat.  If
% a string S is given, loads from S.mat.  You can pass the
% trailing .mat as well and it will still work.
%
% If whichCal is specified, the whichCal'th calibration
% in the file is returned.  If whichCal > nCals, an
% empty calibration is returned.  whichCal defaults
% to the most recent calibration.
%
% If the specified file cannot be found, returns empty matrix.
%
% The returned variable cal is a structure containing calibration 
% information.
%
% The returned cell array cals contains all of the calibrations
% stored in the file
%
% The returned string fullFilename is the full path to the calibration
% file.
%
% See also SaveCalFile, CalDataFolder.
%
% 5/28/96  dgp  Wrote it.
% 6/6/96   dgp  Use CalibrationsFolder.
% 6/6/96   dgp  Use whole path in filename so Matlab will only look there.
% 7/25/96  dgp  Use CalDataFolder.
% 8/4/96   dhb  More flexible filename interface.
% 8/21/97  dhb  Rewrite for calibrations stored as cell array.
%               Optional return of entire calibration history.
% 8/26/97  dhb  Handle case of isempty(cals).
%               Added whichCal argument.
% 5/18/99  dhb  Added dir argument.
% 8/15/00  dhb  Modify to handle local/demo cal directories.
% 4/2/13   dhb  Updated for subdir searching logic.
% 6/2/13   dhb  More robust about whether passed filespec contains the trailing '.mat'.

% Get whichCal
if nargin < 2 || isempty(whichCal)
	whichCal = Inf;
end

% Set filespec
if (nargin < 1 || isempty(filespec))
	filename = ['default.mat'];
elseif (ischar(filespec))
    if (~strcmp(filespec(end-3:end),'.mat'))
        filename = [filespec '.mat'];
    else
        filename = filespec;
    end
else
	filename = [sprintf('screen%d.mat', filespec)];
end

% Set the directory
if nargin < 3 || isempty(dir)
	useDir = CalDataFolder(0,filename);
else
	useDir = CalDataFolder(0,filename,dir);
end
fullFilename = [useDir filename];

% If the file doesn't exist in the usual location, take a look in the
% secondary location.
if (~exist(fullFilename, 'file') && (nargin < 3 || isempty(dir)))
	useDir = CalDataFolder(1,filename);
	fullFilename = [useDir filename];
end

% Now read the sucker if it is there.
if exist(fullFilename, 'file')
	eval(['load ' QuoteString(fullFilename)]);
	if isempty(cals) %#ok<NODEF>
		cal = [];
	else
		% Get the number of calibrations.
		nCals = length(cals);
		
		% User the most recent calibration (the last one in the cals cell
		% array) by default.  If the user specified a particular cal file,
		% try to retrieve it or return an empty matrix if the cal index is
		% out of range.
		if whichCal == Inf
			cal = cals{nCals};
		elseif whichCal > nCals || whichCal < 1
			cal = [];
		else
			cal = cals{whichCal};
		end
	end
else
	cal = [];
	cals = {};
end
