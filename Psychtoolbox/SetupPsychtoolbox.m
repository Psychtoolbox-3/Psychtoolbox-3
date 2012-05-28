function SetupPsychtoolbox
% SetupPsychtoolbox - In-place setup of PTB without network access.
%
% This script prepares an already downloaded working copy of Psychtoolbox
% for use with Matlab or Octave. It sets proper paths, performs online
% registration if connected to a network and takes care of special setup
% operations for the Java based GetChar implementation.
%
% This setup routine is meant for people who want to install Psychtoolbox
% but don't have direct access to the internet. Installation in that case
% is a three step procedure:
%
% 1. Download and unpack a full working copy of PTB into your target
% folder. Obviously you need to somehow get a copy, either via conventional
% download from a computer with network connection (See 'help
% DownloadPsychtoolbox' or 'help UpdatePsychtoolbox') or from a helpful
% colleague.
%
% 2. Change your Matlab/Octave working directory to the Psychtoolbox installation
% folder, e.g., 'cd /Applications/Psychtoolbox'.
%
% 3. Type 'SetupPsychtoolbox' to run this script.
%
% Please be aware that the recommended method of installation is via the
% online Subversion system, i.e., DownloadPsychtoolbox and
% UpdatePsychtoolbox. Some functionality may not work with a copy that is
% set up via this script, e.g., PsychtoolboxVersion may provide incomplete
% version information. Convenient upgrades via UpdatePsychtoolbox may be
% impossible. Download size with this method is much higher as well.
%
% If you get stuck, post your question to the forum by email or web:
% web mailto:psychtoolbox@yahoogroups.com
% web http://groups.yahoo.com/group/psychtoolbox/messages/
% web http://groups.yahoo.com/group/psychtoolbox/post/
% Please specify your full name and the version of your operating system,
% MATLAB or Octave and Psychtoolbox.
%

% History:
%
% 01/25/07 mk  Created. Derived from DownloadPsychtoolbox. Basically a
%              stripped down version of it.
%
% 05/31/09 mk  Add support for Octave-3.
% 06/07/09 mk  Apply smallish fix to fix counting, contributed by Chuan Zeng.
% 05/27/12 mk  Strip backwards compatibility support to pre-R2007a.
%              Remove PowerPC support.

% Flush all MEX files: This is needed at least on M$-Windows to
% work if Screen et al. are still loaded.
clear mex

% Check if this is a 64-bit Matlab or Octave, which we don't support at all:
if strcmp(computer,'PCWIN64') || ...
  (~isempty(strfind(computer, '_64')) && isempty(strfind(computer, 'linux')) && isempty(strfind(computer, 'apple'))) 
    fprintf('Psychtoolbox does not work on a 64 bit version of Matlab or Octave on MS-Windows.\n');
    fprintf('You need to install a 32 bit Matlab or Octave to install and use Psychtoolbox.\n');
    fprintf('Use with 64 bit Matlab or (soon) Octave is fully supported on GNU/Linux and MacOSX.\n');
    fprintf('ERROR: See also http://psychtoolbox.org/wikka.php?wakka=Faq64BitSupport.\n');
    error('Tried to setup on a 64 bit version of Matlab or Octave, which is not supported on this operating system.');
end

if strcmp(computer,'MAC')
    fprintf('This version of Psychtoolbox is no longer supported under MacOSX on the Apple PowerPC hardware platform.\n');
    fprintf('You can get modern versions of Psychtoolbox-3 for Linux if you choose to install GNU/Linux on your PowerPC\n');
    fprintf('machine. These are available from the GNU/Debian project and a future Ubuntu 12.10 release\n.');
    fprintf('Alternatively you can download old - totally unsupported - releases of Psychtoolbox version 3.0.9 from GoogleCode:\n');
    fprintf('http://code.google.com/p/psychtoolbox-3/ \n\n');
    error('Apple MacOSX on Apple PowerPC computers is no longer supported by this Psychtoolbox version.');
end

% Check OS
isWin = strcmp(computer,'PCWIN') || strcmp(computer, 'i686-pc-mingw32');
isOSX = ~isempty(strfind(computer, 'MAC')) || ~isempty(strfind(computer, 'apple-darwin'));
isLinux = strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64') || ~isempty(strfind(computer, 'linux-gnu'));

if ~isWin && ~isOSX && ~isLinux
    os = computer;
    if strcmp(os,'MAC2')
        os = 'Mac OS9';
    end
    fprintf('Sorry, this updater doesn''t support your operating system: %s.\n', os);
    fprintf([mfilename ' can only install the new (OSX, Linux and Windows) \n'...
        'OpenGL-based versions of the Psychtoolbox-3. To install the older (OS9 and Windows) \n'...
        'versions (not based on OpenGL, aka PTB-2) please go to the legacy Psychtoolbox website: \n'...
        'web http://psychtoolbox.org/PTB-2/index.html\n']);
    error(['Your operating system is not supported by ' mfilename '.']);
end

% Locate ourselves:
targetdirectory=fileparts(mfilename('fullpath'));
if ~strcmpi(targetdirectory, pwd)
    error('You need to change your working directory to the Psychtoolbox folder before running this routine!');
end

fprintf('Will setup working copy of the Psychtoolbox folder inside: %s\n',targetdirectory);
fprintf('\n');

% Does SAVEPATH work?
if exist('savepath') %#ok<EXIST>
   err=savepath;
else
   err=path2rc;
end

if err
    try
        % If this works then we're likely on Matlab:
        p=fullfile(matlabroot,'toolbox','local','pathdef.m');
        fprintf(['Sorry, SAVEPATH failed. Probably the pathdef.m file lacks write permission. \n'...
            'Please ask a user with administrator privileges to enable \n'...
            'write by everyone for the file:\n\n''%s''\n\n'],p);
    catch %#ok<CTCH>
        % Probably on Octave:
        fprintf(['Sorry, SAVEPATH failed. Probably your ~/.octaverc file lacks write permission. \n'...
            'Please ask a user with administrator privileges to enable \n'...
            'write by everyone for that file.\n\n']);
    end
    
    fprintf(['Once "savepath" works (no error message), run ' mfilename ' again.\n']);
    fprintf('Alternatively you can choose to continue with installation, but then you will have\n');
    fprintf('to resolve this permission isssue later and add the path to the Psychtoolbox manually.\n\n');
    answer=input('Do you want to continue the installation despite the failure of SAVEPATH (yes or no)? ','s');
    if ~strcmp(answer,'yes')
        fprintf('\n\n');
        error('SAVEPATH failed. Please get an administrator to allow everyone to write pathdef.m.');
    end
end

% Handle Windows ambiguity of \ symbol being the filesep'arator and a
% parameter marker:
if isWin
    searchpattern = [filesep filesep 'Psychtoolbox[' filesep pathsep ']'];
    searchpattern2 = [filesep filesep 'Psychtoolbox'];
else
    searchpattern  = [filesep 'Psychtoolbox[' filesep pathsep ']'];
    searchpattern2 = [filesep 'Psychtoolbox'];
end

% Remove "Psychtoolbox" from path:
while any(regexp(path, searchpattern))
    fprintf('Your old Psychtoolbox appears in the MATLAB/OCTAVE path:\n');
    paths=regexp(path,['[^' pathsep ']+'],'match');
    fprintf('Your old Psychtoolbox appears %d times in the MATLAB/OCTAVE path.\n',length(paths));
    % Old and wrong, counts too many instances: fprintf('Your old Psychtoolbox appears %d times in the MATLAB/OCTAVE path.\n',length(paths));
    answer=input('Before you decide to delete the paths, do you want to see them (yes or no)? ','s');
    if ~strcmp(answer,'yes')
        fprintf('You didn''t say "yes", so I''m taking it as no.\n');
    else
        for p=paths
            s=char(p);
            if any(regexp(s,searchpattern2))
                fprintf('%s\n',s);
            end
        end
    end
    answer=input('Shall I delete all those instances from the MATLAB/OCTAVE path (yes or no)? ','s');
    if ~strcmp(answer,'yes')
        fprintf('You didn''t say yes, so I cannot proceed.\n');
        fprintf('Please use the MATLAB "File:Set Path" command or its Octave equivalent to remove all instances of "Psychtoolbox" from the path.\n');
        error('Please remove Psychtoolbox from MATLAB/OCTAVE path.');
    end
    for p=paths
        s=char(p);
        if any(regexp(s,searchpattern2))
            % fprintf('rmpath(''%s'')\n',s);
            rmpath(s);
        end
    end
    if exist('savepath') %#ok<EXIST>
       savepath;
    else
       path2rc;
    end

    fprintf('Success.\n\n');
end

% Add Psychtoolbox to MATLAB/OCTAVE path
fprintf('Now adding the new Psychtoolbox folder (and all its subfolders) to your MATLAB/OCTAVE path.\n');
p=targetdirectory;
pp=genpath(p);
addpath(pp);

if exist('savepath') %#ok<EXIST>
   err=savepath;
else
   err=path2rc;
end

if err
    fprintf('SAVEPATH failed. Psychtoolbox is now already installed and configured for use on your Computer,\n');
    fprintf('but i could not save the updated Matlab/Octave path, probably due to insufficient permissions.\n');
    fprintf('You will either need to fix this manually via use of the path-browser (Menu: File -> Set Path),\n');
    fprintf('or by manual invocation of the savepath command (See help savepath). The third option is, of course,\n');
    fprintf('to add the path to the Psychtoolbox folder and all of its subfolders whenever you restart Matlab.\n\n\n');
else 
    fprintf('Success.\n\n');
end

fprintf(['Now setting permissions to allow everyone to write to the Psychtoolbox folder. This will \n'...
    'allow future updates by every user on this machine without requiring administrator privileges.\n']);

try
    if isOSX || isLinux
        [s,m]=fileattrib(p,'+w','a','s'); % recursively add write privileges for all users.
    else
        [s,m]=fileattrib(p,'+w','','s'); % recursively add write privileges for all users.
    end
catch %#ok<CTCH>
    s = 0;
    m = 'Setting file attributes is not supported under Octave.';
end

if s
    fprintf('Success.\n\n');
else
    fprintf('\nFILEATTRIB failed. Psychtoolbox will still work properly for you and other users, but only you\n');
    fprintf('or the system administrator will be able to run the UpdatePsychtoolbox script to update Psychtoolbox,\n');
    fprintf('unless you or the system administrator manually set proper write permissions on the Psychtoolbox folder.\n');
    fprintf('The error message of FILEATTRIB was: %s\n\n', m);
end

if exist('PsychtoolboxPostInstallRoutine.m', 'file')
   % Notify the post-install routine of the "pseudo-update" It will
   % determine the proper flavor by itself.
   PsychtoolboxPostInstallRoutine(1);
end

% Puuh, we are done :)
return
