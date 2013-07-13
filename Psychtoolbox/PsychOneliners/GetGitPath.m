function gitpath = GetGitPath
% gitpath = GetGitPath -- Return auto-detected installation path
% for git client, if any. Return empty string if auto-detection not
% possible. Typical usage is like this:
%
% mygitcommand = [GetGitPath 'git describe']; system(mygitcommand);
%
% GetGitPath will return the path to be prefixed in front of the git
% executable. If none can be found, the git executable will be executed
% without path spec. If it is installed in the system executable search
% path, it will then still work.
%
% The function simply checks if the git executable is in the Matlab path
% and returns a proper path-spec. If it isn't found in the Matlab path, it
% tries default path locations for OS-X and Windows. If that doesn't work,
% it returns an empty string.

% History:
% 07/11/13 Written, based on GetSubversionPath (DHB).

% Check for alternative install location of Git:
if IsWin
	% Search for Windows executable in Matlab's path:
	gitpath = which('git.exe');
else
	% Search for Unix executable in Matlab's path:
	gitpath = which('git.');
end

% Found one?
if ~isempty(gitpath)
	% Extract basepath and use it:
	gitpath=[fileparts(gitpath) filesep];
else
	% Could not find git executable in Matlabs path. Check the default
	% install location on OS-X and abort if it isn't there. On M$-Win we
	% simply have to hope that it is in some system dependent search path.

	% Currently, we only know how to check this for Mac OSX.
	if IsOSX
		gitpath = '';
        
		if isempty(gitpath) && exist('/usr/bin/git','file')
			gitpath='/usr/bin/';
		end

		if isempty(gitpath) && exist('/usr/local/bin/git','file')
			gitpath='/usr/local/bin/';
		end

		if isempty(gitpath) && exist('/bin/git','file')
			gitpath='/bin/';
		end

		if isempty(gitpath) && exist('/opt/local/bin/git', 'file')
			gitpath = '/opt/local/bin/';
		end
	end
end

return;
