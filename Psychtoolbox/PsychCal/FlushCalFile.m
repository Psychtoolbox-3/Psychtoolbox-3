function FlushCalFile(filespec,nKeep)
% FlushCalFile([filespec],[nKeep])
%
% Flush all but the most recent calibrations in a file.
%
% If no filespec is given, flushes file default.mat.  If
% a an integer is given, flushes file screenN.mat.  If
% a string is given, loads from string.mat.
%
% Argument nKeep specifies how many calibrations to
% keep.  If fewer than nKeep are in the file, all are
% kept.  If nKeep is zero, file is emptied.  nKeep
% defaults to 1.
%
% If specified file cannot be found does nothing.
%
% 8/26/97  dhb  Wrote it.
% 8/21/00  dhb  Update for dual cal dir scheme.  Not tested hard.

% Set nKeep
if (nargin < 2 || isempty(nKeep))
	nKeep = 1;
end

% Set the filename
if (nargin < 1 || isempty(filespec))
	filename = [CalDataFolder 'default.mat'];
elseif (ischar(filespec))
	filename = [CalDataFolder filespec '.mat'];
else
	filename = [CalDataFolder sprintf('screen%d.mat',filespec)];
end

% Make sure file is present before calling load
file=fopen(filename);

% If not, make sure to try secondary directory.
if (file == -1 && (nargin < 3 || isempty(dir)))
	useDir = CalDataFolder(1);
	if (nargin < 1 || isempty(filespec))
		filename = [useDir 'default.mat'];
	elseif (ischar(filespec))
		filename = [useDir filespec '.mat'];
	else
		filename = [useDir sprintf('screen%d.mat',filespec)];
	end
	file = fopen(filename);
end

% Now read the sucker if it is there.
if (file ~= -1)
	fclose(file);
	eval(['load ' QuoteString(filename)])
	nCals = size(cals,2);
	cal = cals{nCals};
	if (nCals < nKeep)
		return;
	elseif (nKeep == 0)
		cals = {};
	elseif (nKeep == 1)
		cals = {cals{nCals}};
	else
		cals = {cals{nCals-nKeep+1:nCals}};
	end

	% Save the flushed calibration file
	eval(['save ' QuoteString(filename) ' cals']);
end
