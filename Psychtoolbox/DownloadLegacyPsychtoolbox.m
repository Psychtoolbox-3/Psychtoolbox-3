function DownloadLegacyPsychtoolbox(targetdirectory,downloadmethod,targetRevision,flavor)
% DownloadLegacyPsychtoolbox([targetdirectory] [,downloadmethod=0] [,targetRevision][,flavor])
%
% CAUTION: This script is for legacy downloads. Use DownloadPsychtoolbox
% for downloads of the current Psychtoolbox 3.0.10 or later. Psychtoolbox
% 3.0.9 and earlier are completely unsupported and unmaintained!
%
% This script downloads *old versions* of Psychtoolbox-3, specifically
% version 3.0.9 or earlier from the Subversion master server at GoogleCode
% to your disk, creating your working copy, ready to use as a new toolbox
% in your MATLAB/OCTAVE application. Subject to your permission, any old
% installation of the Psychtoolbox is first removed. It's a careful
% program, checking for all required resources and privileges before it
% starts.
%
% CAUTION: Psychtoolbox 3.0.9 and earlier *will not work* with 64 bit
% versions of Matlab or Octave, except if you use a GNU/Linux system, e.g.,
% Debian GNU/Linux or Ubuntu Linux 10.04 or later. The NeuroDebian project
% (see http://neuro.debian.net) provides a very convinient installation of
% Psychtoolbox for both 32 bit and 64 bit versions of Octave and Matlab via
% the regular package management system of your Linux distribution.
%
% On Mac OSX, all parameters are optional. On MS-Windows and GNU/Linux, the
% first parameter "targetdirectory" with the path to the installation
% target directory is required. The "targetdirectory" name may not contain
% any white space, otherwise download will fail with mysterious error
% messages!
% 
% On OSX, your working copy of the Psychtoolbox will be placed in either
% your /Applications or your /Users/Shared folder (depending on permissions
% and your preference), or you may specify a 'targetdirectory', as you
% prefer.
%
% On Microsoft Windows, you must specify the full path, including
% the drive name where Psychtoolbox should be installed, e.g.,
% 'C:\MyToolboxes\'.
%
% The desired flavor of a Psychtoolbox release can be selected via the
% optional "flavor" parameter: By default, 'beta' (aka 'current') will be
% installed if you don't specify otherwise, as this is almost always the
% best possible choice. You can download an old versioned release via a
% namestring like 'Psychtoolbox-x.y.z', e.g., 'Psychtoolbox-3.0.7' if you'd
% want to download version 3.0.7. This is only useful if you run a very old
% operating system or Matlab version that isn't supported by the current
% "beta" anymore, so you'd need to stick with an old versioned release.
% People that really love trouble can also download the 'unsupported'
% flavor. "Unsupported" was formerly known as "stable" but the new name
% reflects reality much better and accurately describes the level of support
% you can expect from us if you use it and run into any trouble!
%
%
% Normally your download should just work(TM). The installer knows three
% different methods of download and tries all of them if neccessary, ie.,
% if the preferred method fails, the 2nd best is tried etc. Should the
% installer get stuck for an inappropriate amount of time (More than 5-10
% minutes), you can try to abort it and restart it, providing the
% additional 'downloadmethod' parameter with a setting of either 0 or 1,
% to change the order of tried download methods to prevent the downloader
% from getting stuck with a specific method in rare cases. Very
% infrequently, the download servers may be overloaded or down for
% maintenance, resulting in download failure. In that case, please retry a
% few hours later.
%
%
% The "targetRevision" argument is optional and should be normally omitted.
% Normal behaviour is to download the latest revision of Psychtoolbox.
% If you provide a specific targetRevision, then this script will
% install a copy of Psychtoolbox according to the specified revision.
%
% This is only useful if you experience problems and want
% to revert to an earlier known-to-be-good release.
%
% Revisions can be specified by a revision number, a specific date, or by
% the special flag 'PREV' which will choose the revision before the
% most current one.
%
%
% INSTALLATION INSTRUCTIONS: The Wiki contains much more up to date
% instructions. If in doubt, follow instructions on the Wiki!
%
% 1. If you don't already have it, you must install the Subversion client.
% For Mac OSX, download the latest Mac OSX Subversion client from:
% web http://metissian.com/projects/macosx/subversion/
% (You can ignore the Subversion README file. If you do read it, you can
% skip the instruction to manually add /usr/local/bin to your unix path.
% That's tricky to do, and not needed for installation and updates because
% we always specify the full path.) Please note that OS/X 10.5 "Leopard"
% and later often come already with Subversion preinstalled, so you may
% be able to skip step 1.
%
% For Windows, download the Windows Subversion client from:
% web http://subversion.tigris.org/files/documents/15/34093/svn-1.4.0-setup.exe
% Install the Subversion client on your machine by double-clicking the
% installer and following the instructions. After installation of the
% Subversion client, you will need to exit and restart Matlab or Octave, so it
% can find the new subversion executable. In many cases it may be
% neccessary to even reboot your computer after installation of subversion.
% Btw. you should avoid to install the client into a path that contains
% blanks/spaces/white-space as this can lead to download failures in some
% cases, e.g., 'C:\Program Files\...' may be bad because there is a blank
% between the "Program" and "Files".
%
% Alternatively, if you don't have the neccessary permissions to install
% Subversion into a system folder, you can install Subversion into an
% arbitrary folder on your system (excluding ones with blanks in their
% path) and then add that folder to your Matlab or Octave path. E.g. you installed
% into D:\MyOwnFolder\Subversion\ . Then you can do this:
% addpath('D:\MyOwnFolder\Subversion\'). Our installer should find the
% client then.
%
% 2. On MacOS/X, to install the Psychtoolbox in the default location
% (/Applications or, failing that, /Users/Shared). Just type:
%
% DownloadPsychtoolbox
%
% Our standard option is in the Applications folder, but note that, as with
% installation of any software, you'll need administrator privileges. Also
% note that if you put the toolbox in the Applications folder, you'll need
% to reinstall it when MATLAB / OCTAVE is updated on your machine. If you must
% install without access to an administrator, we offer the option of
% installing into the /Users/Shared/ folder instead. If you must install
% the Psychtoolbox in some other folder, then specify it in the optional
% first argument of your call.
%
% On Windows or Linux, provide a pathname, e.g.:
% DownloadPsychtoolbox('C:\MyToolboxes\');
%
% That's it. Any pre-existing installation of the Psychtoolbox will be
% removed (if you approve). The program will then download the latest
% Psychtoolbox and update your MATLAB / OCTAVE path and other relevant system settings.
%
% Enjoy! If you're new to this, you might start by typing "help
% Psychtoolbox".
%
% P.S. If you get stuck, first check the FAQ section and Download section of
% our Wiki at http://www.psychtoolbox.org. If that doesn't help, post your
% question to the forum by email or web:
%
% web mailto:psychtoolbox@yahoogroups.com
% web http://groups.yahoo.com/group/psychtoolbox/messages/
% web http://groups.yahoo.com/group/psychtoolbox/post/
%
% Please specify your full name and the version of your operating system,
% MATLAB / OCTAVE, and psychtoolbox.
%
% UPGRADE INSTRUCTIONS:
%
% To upgrade your copy of Psychtoolbox, at any time, to incorporate the
% latest bug fixes, enhancements, and new features, just type:
% UpdatePsychtoolbox
% 
% UpdatePsychtoolbox cannot change the flavor of your
% Psychtoolbox. To change the flavor, run DownloadPsychtoolbox to
% completely discard your old installation and get a fresh copy with the
% requested flavor.
% 
% PERMISSIONS:
%
% There's a thorny issue with permissions on OS/X. It may not be possible to
% install into /Applications (or whatever the targetdirectory is) with the
% user's existing privileges. The normal situation on Mac OSX is that a few
% users have "administrator" privileges, and many don't. By default,
% writing to the /Applications folder requires administrator privileges.
%
% Thus all OSX installers routinely demand an extra authorization (if
% needed), asking the user to type in the name and password of an
% administrator before proceeding. We haven't yet figured out how to do
% that, but we want to offer that option. This conforms to normal
% installation of an application under Mac OS X.
%
% DownloadPsychtoolbox creates the Psychtoolbox folder with permissions set
% to allow writing by everyone. Our hope is that this will allow updating
% (by UpdatePsychtoolbox) without need for administrator privileges.
%
% Some labs that may want to be able to install without access to an
% administrator. For them we offer the fall back of installing Psychtoolbox
% in /Users/Shared/, instead of /Applications/, because, by default,
% /Users/Shared/ is writeable by all users.
%
% SAVEPATH
%
% Normally all users of MATLAB / OCTAVE use the same path. This path is normally
% saved in MATLABROOT/toolbox/local/pathdef.m, where "MATLABROOT" stands
% for the result returned by running that function in MATLAB, e.g.
% '/Applications/MATLAB.app/Contents/Matlab14.1'. Since pathdef.m is inside
% the MATLAB package, which is normally in the Applications folder,
% ordinary users (not administrators) cannot write to pathdef.m. They'll
% get an error message whenever they try to save the path, e.g. by typing
% "savepath". Most users will find this an unacceptable limitation. The
% solution is very simple, ask an administrator to use File Get Info to set
% the pathdef.m file permissions to allow write by everyone. This needs to
% be done only once, after installing MATLAB.
% web http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_env/ws_pat18.html
%
%

% History:
%
% 11/02/05 mk  Created.
% 11/25/05 mk  Bug fix for 'targetdirectory' provided by David Fencsik.
% 01/13/06 mk  Added support for download of Windows OpenGL-PTB.
% 03/10/06 dgp Expanded the help text, above, incorporating suggestions
%              from Daniel Shima.
% 03/11/06 dgp Check OS. Remove old Psychtoolbox from path and from disk.
%              After downloading, add new Psychtoolbox to path.
% 03/12/06 dgp Polished error message regarding missing svn.
% 03/13/06 dgp Changed default targetdirectory from PWD to ~/Documents/.
% 03/14/06 dgp Changed default targetdirectory to /Applications/, and
%              if not sufficiently privileged, then /Users/Shared/.
%              Check privilege to create folder. Check SAVEPATH.
% 06/05/06 mk  On Windows, we require the user to pass the full path to
%              the installation folder, because there is no well-defined
%              equivalent to the Mac OS X /Applications/ folder. Also, the order
%              of operations was changed to maximize the chance of getting a
%              working PTB installation despite minor failures of commands
%              like savepath, or fileattrib. We allow the user to decide
%              whether to delete her old Psychtoolbox folders or to
%              retain multiple copies of Psychtoolbox (e.g., beta and stable)
%              so each user can choose between Beta and Stable.
%              We no longer copy UpdatePsychtoolbox.m, since it's included in the 
%              new Psychtoolbox folder 
% 06/06/06 dgp Cosmetic editing of comments, above.
% 06/27/06 dgp Cosmetic editing of comments and messages. Check for spaces
%              in targetdirectory name.
% 9/23/06  mk  Add clear mex call to flush mex files before downloading.
% 10/5/06  mk  Add detection code for MacOS-X on Intel Macs.
% 10/28/06 dhb Allow 'current' as a synonym for 'beta'.
% 11/21/06 mk  Allow alternate install location for svn client: Installer
%              will find the svn executable if its installation folder is
%              included in the Matlab path.
% 02/17/07 mk  Convert flavor spec to lower case in case it isn't.
% 03/15/07 mk  Detection code for Windows 64 bit added.
% 07/18/07 mk  Changed default for flavor from 'stable' to 'beta'.
% 09/27/07 mk  Add another fallback path: Download via https protocol to
%              maybe bypass proxy-servers.
% 10/29/07 mk  Small fix for Kerstin Preuschoffs bugreport: Download of old
%              versions didn't work anymore, becuase 'flavor' string was
%              lower-cased.
% 11/17/07 mk  Prioritized https:// checkout over http:// checkout --> More
%              likely to bypass proxy servers and succeed.
% 12/25/07 mk  Add check for white spaces in path to subversion client - Output warning, if so.
% 12/25/07 mk  Add optional 'downloadmethod' parameter, which allows to
%              select initial choie of download protocol to use, in order to allow to
%              bypass misconfigured proxies and firewalls. Problem reported
%              by Patrick Mineault.
% 01/08/08 mk  On OS/X, add an additional search path for 'svn': /usr/bin/
%              as per suggestion of Donald Kalar. Presumably, Apples Leopard ships
%              with a svn client preinstalled in that location.
%
% 05/07/08 mk  Add better handling of default values. Add 'targetRevision'
%              parameter as option, just as in UpdatePsychtoolbox.m
%
% 06/16/08 mk  Change default initial download protcol from svn: to https:,
%              as Berlios seems to have disabled svnserve protocol :-(
%
% 10/01/08 mk  Add interactive output/query for svn client on the Unices.
%              This to work-around questions of the client about accepting
%              security certificates...
% 01/05/09 mk  Remove && or || to make old Matlab versions happier.
% 01/05/09 mk  Remove && or || to make old Matlab versions happier.
% 01/05/09 mk  Change order of call arguments to make choice of default
%              'beta' flavor more convenient and choice of 'stable'
%              inconvenient. Also add some additional warning text and
%              dialog to tell user how bad 'stable' is and to give him a
%              chance to reconsider.
% 03/22/09 mk  Update help text again. Rename 'stable' into 'unsupported'.
% 05/31/09 mk  Add support for Octave-3.
% 10/05/09 mk  Strip trailing fileseperator from targetDirectory, as
%              suggested by Erik Flister to avoid trouble with new svn
%              clients.
% 12/27/10 mk  Redirect 'beta' downloads on Matlab versions < 6.5 to the
%              special "Psychtoolbox-3.0.8-PreMatlab6.5" compatibility
%              version - The last one to support pre 6.5 Matlab's.
%
% 10/31/11 mk  Change location of SVN repository to
%              http://psychtoolbox-3.googlecode.com/svn/
%              our new home, now that Berlios is shutting down.
%
% 05/30/12 mk  Fix disasterous bug: Answering anything but 'y' to the
%              question if old PTB folder should be deleted, caused
%              *deletion* of the folder! Oh dear! This bug present since
%              late 2011.

% Flush all MEX files: This is needed at least on M$-Windows for SVN to
% work if Screen et al. are still loaded.
clear mex

% Check if this is a 64-bit Matlab or Octave, which we don't support at all:
if strcmp(computer,'PCWIN64') | strcmp(computer,'MACI64') | ...
  (~isempty(findstr(computer, '_64')) & isempty(findstr(computer, 'linux'))) %#ok<OR2>
    fprintf('Psychtoolbox does not work on a 64 bit version of Matlab or Octave on MS-Windows or Apple-OSX.\n');
    fprintf('You need to install a 32 bit Matlab or Octave to install and use Psychtoolbox.\n');
    fprintf('Use with 64 bit Matlab or Octave is fully supported on GNU/Linux.\n');
    fprintf('ERROR: See also http://psychtoolbox.org/wikka.php?wakka=Faq64BitSupport.\n');
    error('Tried to install on a 64 bit version of Matlab or Octave, which is not supported on this operating system.');
end

% Check OS
IsWin=strcmp(computer,'PCWIN') | strcmp(computer,'PCWIN64')| strcmp(computer, 'i686-pc-mingw32');
IsOSX=strcmp(computer,'MAC') | strcmp(computer,'MACI') | ~isempty(findstr(computer, 'apple-darwin'));
IsLinux=strcmp(computer,'GLNX86') | strcmp(computer,'GLNXA64') | ~isempty(findstr(computer, 'linux-gnu'));

if ~IsWin & ~IsOSX & ~IsLinux
    os=computer;
    if strcmp(os,'MAC2')
        os='Mac OS9';
    end
    fprintf('Sorry, this installer doesn''t support your operating system: %s.\n',os);
    fprintf([mfilename ' can only install the new (OSX, Linux and Windows) \n'...
        'OpenGL-based versions of the Psychtoolbox. To install the older (OS9 and Windows) \n'...
        'versions (not based on OpenGL) please go to the psychtoolbox website: \n'...
        'web http://psychtoolbox.org/download.html\n']);
    error(['Your operating system is not supported by ' mfilename '.']);
end

if nargin < 1
    targetdirectory = [];
end

if isempty(targetdirectory)
    if IsOSX
        % Set default path for OSX install:
        targetdirectory=fullfile(filesep,'Applications');
    else
        % We do not have a default path on Windows, so the user must provide it:
        fprintf('You did not provide the full path to the directory where Psychtoolbox should be\n');
        fprintf('installed. This is required for Microsoft Windows and Linux installation. Please enter a full\n');
        fprintf('path as the first argument to this script, e.g. DownloadPsychtoolbox(''C:\\Toolboxes\\'').\n');
        error('For Windows and Linux, the call to %s must specify a full path for the location of installation.',mfilename);
    end     
end

% Strip trailing fileseperator, if any:
if targetdirectory(end) == filesep
    targetdirectory = targetdirectory(1:end-1);
end

% Override for download method provided?
if nargin < 2
    downloadmethod = [];
end

if isempty(downloadmethod)
    % Try 0 by default (http://):
    downloadmethod = 0;
else
    if downloadmethod < 0 | downloadmethod > 1
        error('Invalid downloadmethod provided. Valid are values 0 and 1.');
    end
end

if nargin < 3
    targetRevision = [];
end

if isempty(targetRevision)
    targetRevision = '';
else
    fprintf('Target revision: %s \n', targetRevision);
    targetRevision = [' -r ' targetRevision ' '];
end

% Set flavor defaults and synonyms
if nargin < 4
    flavor = [];
end

if isempty(flavor)
    flavor='beta';
end

% Make sure that flavor is lower-case, unless its a 'Psychtoolbox-x.y.z'
% spec string which is longer than 10 characters and mixed case:
if length(flavor) < 10
    % One of the short flavor spec strings: lowercase'em:
    flavor = lower(flavor);
end

switch (flavor)
    % 'current' is a synonym for 'beta'.
    case 'beta'
    case 'current'
        flavor = 'beta';        
    case 'stable'
        fprintf('\n\n\nYou request download of the "stable" flavor of Psychtoolbox.\n');
        fprintf('The "stable" flavor is no longer available, it has been renamed to "unsupported".\n');
        fprintf('If you really want to use the former "stable" flavor, please retry the download\n');
        fprintf('under the new name "unsupported".\n\n');
        error('Flavor "stable" requested. This is no longer available.');
    case 'unsupported'
        % Very bad choice! Give user a chance to reconsider...
        fprintf('\n\n\nYou request download of the "unsupported" flavor of Psychtoolbox.\n');
        fprintf('Use of the "unsupported" flavor is strongly discouraged! It is outdated and contains\n');
        fprintf('many bugs and deficiencies that have been fixed in the recommended "beta" flavor years ago.\n');
        fprintf('"unsupported" is no longer maintained and you will not get any support if you have problems with it.\n');
        fprintf('Please choose "beta" unless you have very good reasons not to do so.\n\n');
        fprintf('If you answer "no" to the following question, i will download the recommended "beta" flavor instead.\n');
        answer=input('Do you want to continue download of "unsupported" flavor despite the warnings (yes or no)? ','s');
        if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
            flavor = 'beta';
            fprintf('Download of "unsupported" flavor cancelled, will download recommended "beta" flavor instead...\n');
        else
            fprintf('Download of "unsupported" flavor proceeds. You are in for quite a bit of pain...\n');            
        end

        fprintf('\n\nPress any key to continue...\n');
        pause;
        
    otherwise
        fprintf('\n\n\nHmm, requested flavor is the unusual flavor: %s\n',flavor);
        fprintf('Either you request something exotic, or you made a typo?\n');
        fprintf('We will see. If you get an error, this might be the first thing to check.\n');
        fprintf('Press any key to continue...\n');
        pause;
end

if strcmp(flavor, 'beta')
    % Check if this is Matlab of version prior to V 6.5:
    v = ver('matlab');
    if ~isempty(v)
        v = v(1).Version; v = sscanf(v, '%i.%i.%i');
        if (v(1) < 6) | ((v(1) == 6) & (v(2) < 5)) %#ok<AND2,OR2>
            % Matlab version < 6.5 detected. This is no longer
            % supported by current PTB beta. Redirect to the last
            % functional PTB for such ancient Matlab's:
            flavor = 'Psychtoolbox-3.0.8-PreMatlab6.5';
            fprintf('\n\n\nYou request download of the "beta" flavor of Psychtoolbox.\n');
            fprintf('The "beta" flavor is no longer available for your version of Matlab.\n');
            fprintf('Current "beta" only works on Matlab Version 6.5 or later.\n');
            fprintf('I will download "%s" instead, which is the last working\n', flavor);
            fprintf('version of Psychtoolbox for older Matlab releases. Please consider\n');
            fprintf('upgrading to a recent Matlab version or switching to GNU/Octave 3.2.x. Both\n');
            fprintf('should provide better support, performance and a richer feature set.\n\n');
            fprintf('Running on your ancient Matlab should work, but is not supported\n');
            fprintf('anymore. If you run into any problems or bugs, you are on your own.\n');

            fprintf('\n\nPress any key to continue installation or abort now via CTRL+C...\n\n\n');
            pause;
        end
    end
end

fprintf('DownloadPsychtoolbox(''%s'',''%s'')\n',targetdirectory, flavor);
fprintf('Requested flavor is: %s\n',flavor);
fprintf('Requested location for the Psychtoolbox folder is inside: %s\n',targetdirectory);
fprintf('\n');

% Check for alternative install location of Subversion:
if IsWin
    % Search for Windows executable in path:
    svnpath = which('svn.exe');
else
    % Search for Unix executable in path:
    svnpath = which('svn.');
end

% Found one?
if ~isempty(svnpath)
    % Extract basepath and use it:
    svnpath=[fileparts(svnpath) filesep];
else
    % Could not find svn executable in path. Check the default
    % install location on OS-X and abort if it isn't there. On M$-Win we
    % simply have to hope that it is in some system dependent search path.

    % Currently, we only know how to check this for Mac OSX.
    if IsOSX
        % Try OS/X 10.5 Leopard install location for svn first:
        svnpath='/usr/bin/';
        if exist('/usr/bin/svn','file')~=2
            % This would have been the default install location of the svn
            % client bundled with OS/X 10.5.x Leopard. Let's try the
            % default install location from the web installer:
            svnpath='/usr/local/bin/';
            if exist('/usr/local/bin/svn','file')~=2
                fprintf('The Subversion client "svn" is not in its expected\n');
                fprintf('location "/usr/local/bin/svn" on your disk. Please \n');
                fprintf('download and install the most recent Subversion client from:\n');
                fprintf('web http://metissian.com/projects/macosx/subversion/ -browser\n');
                fprintf('and then run %s again.\n',mfilename);
                error('Subversion client is missing. Please install it.');
            end
        end
    end
end

if ~isempty(svnpath)
    fprintf('Will use the svn client which is located in this folder: %s\n', svnpath);
end

if any(isspace(svnpath))
    fprintf('WARNING! There are spaces (blanks) in the path to the svn client executable (see above).\n');
    fprintf('On some systems this can cause a download failure, with some error message that may look\n');
    fprintf('roughly like this: %s is not recognized as an internal or external command,\n', svnpath(1:min(find(isspace(svnpath)))));
    fprintf('operable program or batch file.\n\n');
    fprintf('Should the download fail with such a message then move/install the svn.exe program into a\n');
    fprintf('folder whose path does not contain any blanks/spaces and retry.\n\n');
    warning('Spaces in path to subversion client -- May cause download failure.');
end

% Does SAVEPATH work?
if exist('savepath')
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
    catch
        % Probably on Octave:
        fprintf(['Sorry, SAVEPATH failed. Probably your ~/.octaverc file lacks write permission. \n'...
            'Please ask a user with administrator privileges to enable \n'...
            'write by everyone for that file.\n\n']);
    end
    
    fprintf(['Once "savepath" works (no error message), run ' mfilename ' again.\n']);
    fprintf('Alternatively you can choose to continue with installation, but then you will have\n');
    fprintf('to resolve this permission isssue later and add the path to the Psychtoolbox manually.\n\n');
    answer=input('Do you want to continue the installation despite the failure of SAVEPATH (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('\n\n');
        error('SAVEPATH failed. Please get an administrator to allow everyone to write pathdef.m.');
    end
end

% Do we have sufficient privileges to install at the requested location?
p='Psychtoolbox123test';
[success,m,mm]=mkdir(targetdirectory,p);
if success
    rmdir(fullfile(targetdirectory,p));
else
	fprintf('Write permission test in folder %s failed.\n', targetdirectory);
    if strcmp(m,'Permission denied')
        if IsOSX
            fprintf([
            'Sorry. You would need administrator privileges to install the \n'...
            'Psychtoolbox into the ''%s'' folder. You can either quit now \n'...
            '(say "no", below) and get a user with administrator privileges to run \n'...
            'DownloadPsychtoolbox for you, or you can install now into the \n'...
            '/Users/Shared/ folder, which doesn''t require special privileges. We \n'...
            'recommend installing into the /Applications/ folder, because it''s the \n'...
            'normal place to store programs. \n\n'],targetdirectory);
            answer=input('Even so, do you want to install the Psychtoolbox into the \n/Users/Shared/ folder (yes or no)? ','s');
            if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
                fprintf('You didn''t say yes, so I cannot proceed.\n');
                error('Need administrator privileges for requested installation into folder: %s.',targetdirectory);
            end
            targetdirectory='/Users/Shared';
        else
            % Windows: We simply fail in this case:
            fprintf([
            'Sorry. You would need administrator privileges to install the \n'...
            'Psychtoolbox into the ''%s'' folder. Please rerun the script, choosing \n'...
            'a location where you have write permission, or ask a user with administrator \n'...
            'privileges to run DownloadPsychtoolbox for you.\n\n'],targetdirectory);
            error('Need administrator privileges for requested installation into folder: %s.',targetdirectory);
        end
    else
        error(mm,m);
    end
end
fprintf('Good. Your privileges suffice for the requested installation into folder %s.\n\n',targetdirectory);

% Delete old Psychtoolbox
skipdelete = 0;
while (exist('Psychtoolbox','dir') | exist(fullfile(targetdirectory,'Psychtoolbox'),'dir')) & (skipdelete == 0)
    fprintf('Hmm. You already have an old Psychtoolbox folder:\n');
    p=fullfile(targetdirectory,'Psychtoolbox');
    if ~exist(p,'dir')
        p=fileparts(which(fullfile('Psychtoolbox','Contents.m')));
        if length(p)==0
            w=what('Psychtoolbox');
            p=w(1).path;
        end
    end
    fprintf('%s\n',p);
    fprintf('That old Psychtoolbox should be removed before we install a new one.\n');
    if ~exist(fullfile(p,'Contents.m'))
        fprintf(['WARNING: Your old Psychtoolbox folder lacks a Contents.m file. \n'...
            'Maybe it contains stuff you want to keep. Here''s a DIR:\n']);
        dir(p)
    end

    fprintf('First we remove all references to "Psychtoolbox" from the MATLAB / OCTAVE path.\n');
    pp=genpath(p);
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(pp);
    warning('on','MATLAB:rmpath:DirNotFound');
    
    if exist('savepath')
       savepath;
    else
       path2rc;
    end

    fprintf('Success.\n');

    answer=input('Shall I delete the old Psychtoolbox folder and all its contents \n(recommended in most cases), (yes or no)? ','s');
    if strcmpi(answer,'yes') || strcmpi(answer,'y')
        skipdelete = 0;
        fprintf('Now we delete "Psychtoolbox" itself.\n');
        [success,m,mm]=rmdir(p,'s');
        if success
            fprintf('Success.\n\n');
        else
            fprintf('Error in RMDIR: %s\n',m);
            fprintf('If you want, you can delete the Psychtoolbox folder manually and rerun this script to recover.\n');
            error(mm,m);
        end
    else
        skipdelete = 1;
    end
end

% Handle Windows ambiguity of \ symbol being the filesep'arator and a
% parameter marker:
if IsWin
    searchpattern = [filesep filesep 'Psychtoolbox[' filesep pathsep ']'];
    searchpattern2 = [filesep filesep 'Psychtoolbox'];
else
    searchpattern  = [filesep 'Psychtoolbox[' filesep pathsep ']'];
    searchpattern2 = [filesep 'Psychtoolbox'];
end

% Remove "Psychtoolbox" from path
while any(regexp(path,searchpattern))
    fprintf('Your old Psychtoolbox appears in the MATLAB / OCTAVE path:\n');
    paths=regexp(path,['[^' pathsep ']*' pathsep],'match');
    fprintf('Your old Psychtoolbox appears %d times in the MATLAB / OCTAVE path.\n',length(paths));
    answer=input('Before you decide to delete the paths, do you want to see them (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('You didn''t say "yes", so I''m taking it as no.\n');
    else
        for p=paths
            s=char(p);
            if any(regexp(s,searchpattern2))
                fprintf('%s\n',s);
            end
        end
    end
    answer=input('Shall I delete all those instances from the MATLAB / OCTAVE path (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('You didn''t say yes, so I cannot proceed.\n');
        fprintf('Please use the MATLAB "File:Set Path" command to remove all instances of "Psychtoolbox" from the path.\n');
        error('Please remove Psychtoolbox from MATLAB / OCTAVE path.');
    end
    for p=paths
        s=char(p);
        if any(regexp(s,searchpattern2))
            % fprintf('rmpath(''%s'')\n',s);
            rmpath(s);
        end
    end
    if exist('savepath')
       savepath;
    else
       path2rc;
    end

    fprintf('Success.\n\n');
end

% Download Psychtoolbox
if IsOSX
    fprintf('I will now download the latest Psychtoolbox for OSX.\n');
else
    if IsLinux
        fprintf('I will now download the latest Psychtoolbox for Linux.\n');
    else
        fprintf('I will now download the latest Psychtoolbox for Windows.\n');
    end
end
fprintf('Requested flavor is: %s\n',flavor);
fprintf('Target folder for installation: %s\n',targetdirectory);
p=fullfile(targetdirectory,'Psychtoolbox');

% Create quoted version of 'p'ath, so blanks in path are handled properly:
pt = strcat('"',p,'"');

% Choose initial download method. Defaults to zero, ie. http protocol:
if downloadmethod < 1
    checkoutcommand=[svnpath 'svn checkout ' targetRevision ' http://psychtoolbox-3.googlecode.com/svn/' flavor '/Psychtoolbox/ ' pt];
else
    % Good to get through many firewalls and proxies:
    checkoutcommand=[svnpath 'svn checkout ' targetRevision ' https://psychtoolbox-3.googlecode.com/svn/' flavor '/Psychtoolbox/ ' pt];
end

fprintf('The following CHECKOUT command asks the Subversion client to \ndownload the Psychtoolbox:\n');
fprintf('%s\n',checkoutcommand);
fprintf('Downloading. It''s nearly 100 MB, which can take many minutes. \nAlas there may be no output to this window to indicate progress until the download is complete. \nPlease be patient ...\n');
fprintf('If you see some message asking something like "accept certificate (p)ermanently, (t)emporarily? etc."\n');
fprintf('then please press the p key on your keyboard, possibly followed by pressing the ENTER key.\n\n');
if IsOSX | IsLinux
    [err]=system(checkoutcommand);
    result = 'For reason, see output above.';
else
    [err,result]=dos(checkoutcommand, '-echo');
end

if err & (downloadmethod == 1)
    % Failed! Let's retry it via http protocol. This may work-around overly
    % restrictive firewalls or otherwise screwed network proxies:
    fprintf('Command "CHECKOUT" failed with error code %d: \n',err);
    fprintf('%s\n\n',result);
    fprintf('Will retry now by use of alternative http protocol...\n');
    checkoutcommand=[svnpath 'svn checkout ' targetRevision ' http://psychtoolbox-3.googlecode.com/svn/' flavor '/Psychtoolbox/ ' pt];
    fprintf('The following alternative CHECKOUT command asks the Subversion client to \ndownload the Psychtoolbox:\n');
    fprintf('%s\n\n',checkoutcommand);
    if IsOSX | IsLinux
        [err]=system(checkoutcommand);
        result = 'For reason, see output above.';
    else
        [err,result]=dos(checkoutcommand, '-echo');
    end    
end

if err & (downloadmethod == 0)
    % Failed! Let's retry it via https protocol. This may work-around overly
    % restrictive firewalls or otherwise screwed network proxies:
    fprintf('Command "CHECKOUT" failed with error code %d: \n',err);
    fprintf('%s\n\n',result);
    fprintf('Will retry now by use of alternative https protocol...\n');
    checkoutcommand=[svnpath 'svn checkout ' targetRevision ' https://psychtoolbox-3.googlecode.com/svn/' flavor '/Psychtoolbox/ ' pt];
    fprintf('The following alternative CHECKOUT command asks the Subversion client to \ndownload the Psychtoolbox:\n');
    fprintf('%s\n\n',checkoutcommand);
    if IsOSX | IsLinux
        [err]=system(checkoutcommand);
        result = 'For reason, see output above.';
    else
        [err,result]=dos(checkoutcommand, '-echo');
    end    
end

if err
    fprintf('Sorry, the download command "CHECKOUT" failed with error code %d: \n',err);
    fprintf('%s\n',result);
    fprintf('The download failure might be due to temporary network or server problems. You may want to try again in a\n');
    fprintf('few minutes. It could also be that the subversion client was not (properly) installed. On Microsoft\n');
    fprintf('Windows you will need to exit and restart Matlab or Octave after installation of the Subversion client. If that\n');
    fprintf('does not help, you will need to reboot your machine before proceeding.\n');
    error('Download failed.');
end
fprintf('Download succeeded!\n\n');

% Add Psychtoolbox to MATLAB / OCTAVE path
fprintf('Now adding the new Psychtoolbox folder (and all its subfolders) to your MATLAB / OCTAVE path.\n');
p=fullfile(targetdirectory,'Psychtoolbox');
pp=genpath(p);
addpath(pp);

if exist('savepath')
   err=savepath;
else
   err=path2rc;
end

if err
    fprintf('SAVEPATH failed. Psychtoolbox is now already installed and configured for use on your Computer,\n');
    fprintf('but i could not save the updated MATLAB / OCTAVE path, probably due to insufficient permissions.\n');
    fprintf('You will either need to fix this manually via use of the path-browser (Menu: File -> Set Path),\n');
    fprintf('or by manual invocation of the savepath command (See help savepath). The third option is, of course,\n');
    fprintf('to add the path to the Psychtoolbox folder and all of its subfolders whenever you restart MATLAB / OCTAVE.\n\n\n');
else 
    fprintf('Success.\n\n');
end

fprintf(['Now setting permissions to allow everyone to write to the Psychtoolbox folder. This will \n'...
    'allow future updates by every user on this machine without requiring administrator privileges.\n']);
try
    if IsOSX | IsLinux
        [s,m]=fileattrib(p,'+w','a','s'); % recursively add write privileges for all users.
    else
        [s,m]=fileattrib(p,'+w','','s'); % recursively add write privileges for all users.
    end
catch
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

fprintf('You can now use your newly installed ''%s''-flavor Psychtoolbox. Enjoy!\n',flavor);
fprintf('Whenever you want to upgrade your Psychtoolbox to the latest ''%s'' version, just\n',flavor);
fprintf('run the UpdatePsychtoolbox script.\n\n');

if exist('PsychtoolboxPostInstallRoutine.m', 'file')
   % Notify the post-install routine of the download and its flavor.
   clear PsychtoolboxPostInstallRoutine;
   PsychtoolboxPostInstallRoutine(0, flavor);
end

% Puuh, we are done :)
return
