function svnInfo = GetSVNInfo(directoryOrFile)
% svnInfo = GetSVNInfo(directoryOrFile)
%
% Description:
% Retrieves the svn information on a specified directory or file.  This is
% essentially a wrapper around the shell command "svn info".
%
% Input:
% directoryOrFile (string) - Directory or file name of interest.
%
% Output:
% svnInfo (struct) - Structure containing the following information:
%   Path
%	URL
%	RepositoryRoot
%	RepositoryUUID
%	Revision
%	NodeKind
%	Schedule
%	LastChangedAuthor
%	LastChangedRev
%	LastChangedDate
%
%	'svnInfo' will be empty if there is no svn info for 'directoryOrFile'.

if nargin ~= 1
	error('Usage: svnInfo = GetSVNInfo(directoryOrFile)');
end

svnInfo = [];

% Look to see if we can find the svn executable on the path.
svnPath = sprintf('%ssvn', GetSubversionPath);
if ~exist(svnPath, 'file')
	fprintf('*** Failed to find svn, returning empty.\n');
	return;
end

% Get the svn info of the specified directory or file.
[status, result] = system(sprintf('%s info "%s"', svnPath, directoryOrFile));
if status ~= 0
	return;
end

if ~isempty(strmatch(result, 'Not a versioned resource'))
	fprintf('*** "%s" is not a versioned resource', directoryOrFile);
	return;
end

% Parse the svn output.
x = textscan(result, '%s', 'Delimiter', '\n');
x = x{1};

for i = 1:length(x)
	if ~isempty(x{i})
		[t, r] = strtok(x{i}, ':');
		svnInfo.(t(t ~= ' ')) = r(3:end);
	end
end
