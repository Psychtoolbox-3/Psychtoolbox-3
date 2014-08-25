function [versionString, versionStructure]=PsychtoolboxVersion
% OS X, Windows, Linux: ___________________________________________________
%
% [versionString, versionStructure]=PsychtoolboxVersion
%
% Return a string identifying this release of the Psychtoolbox.
% The first three numbers identify the base version of Psychtoolbox:
%
% - Leftmost: Increments indicate a disruptive change in the feature
% set and design of the software, by abrupt introduction of design changes.
% This should never happen, as this would mean a completely new product,
% entirely incompatible with the old software.
%
% - Middle: Increments indicate significant enhancements or changes in
% functionality. This will usually only happen every couple of years at
% most.
%
% - Rightmost: A counter to distinguish multiple releases having the same
% leftmost and middle version numbers. This happens if there are backwards
% incompatible changes to the programming interface or functionality which
% may require code adjustments in user code. It also happens if we cancel
% support for platforms (Matlab/Octave versions, operating system versions,
% processor architectures etc.). This happens occassionally.
%
% Numeric values of the three integer fields contained in versionString are
% available in fields of the second return argument, "versionStructure".
%
% The field 'Flavor' defines the subtype of Psychtoolbox being used:
%
% * beta: An experimental release that is already tested by the developers,
% but not yet sufficiently tested or proven in the field. Beta releases
% contain lots of new and experimental features that may be useful to you
% but that may change slightly in behaviour or syntax in the final release,
% making it necessary for you to adapt your code after a software update.
% Beta releases are known to be imperfect and fixing bugs in them is not a
% high priority.  The term 'current' is a synonym for 'beta'. Beta releases
% are the only releases we provide at this point.
%
% * trunk: Development prototypes, for testing and debugging by developers
% and really adventuruous users, not for production use!
%
% * Psychtoolbox-x.y.z: Old, no longer supported Psychtoolbox versions.
%
% * Debian package: A Psychtoolbox provided by GNU/Debian or NeuroDebian.
%
% The revision number and the provided URL allows you to visit the developer
% website in the Internet and get direct access to all development logs
% regarding your working copy of Psychtoolbox.
%
% Be aware that execution of the PsychtoolboxVersion command can take a
% lot of time (in the order of multiple seconds to 1 minute).
%
% Most Psychtoolbox mex files now provide a built-in 'Version' command
% which returns version for themselves.  The version string for the
% built-in version numbers contains a fourth numeric field named "build".
% The build number is a unique serial number.  Mex files distinquished only
% by build numbers were compiled from identical C source files.
%
% _________________________________________________________________________
%
% see also: Screen('Version')

%   2/9/01      awi     added fullfile command for platform-independent pathname
%   6/29/02     dgp     Use new PsychtoolboxRoot function to cope with user-changed folder names.
%   7/12/04     awi     ****** OS X-specific fork from the OS 9 version *******
%                       Noted mex file versioning.
%   7/26/04     awi     Partitioned help and added OS X section.  Enhanced
%                       OS9+Windows section.
%   10/4/05	   awi      Note here that dgp changed "IsWindows" to "IsWin" at unknown date prior
%						between 7/26/04 and 10/4/05.
%
%   5/5/06     mk       Tries to query info from Subversion and displays info about last date
%                       of change, SVN revision and flavor. This code is pretty experimental and
%                       probably also somewhat fragile. And its sloooow.
%   9/17/06    mk       Improvements to parser: We try harder to get the flavor info.
%   10/31/11   mk       Update for our new hoster GoogleCode.
%   04/30/12   mk       Kill MacOS-9 support.
%   05/27/12   mk       Switch over to GitHub hosting.
%   08/06/14   mk       Integrate (Neuro)Debian versioning support. Cleanups.

global Psychtoolbox

if ~isfield(Psychtoolbox,'version')
    Psychtoolbox.version.major=0;
    Psychtoolbox.version.minor=0;
    Psychtoolbox.version.point=0;
    Psychtoolbox.version.string='';
    Psychtoolbox.version.flavor='';
    Psychtoolbox.version.revision=0;
    Psychtoolbox.version.revstring='';
    Psychtoolbox.version.websvn='';

    file=fullfile(PsychtoolboxRoot,'Contents.m');
    f=fopen(file,'r');
    fgetl(f);
    s=fgetl(f);
    fclose(f);
    [cvv,count,errmsg,n]=sscanf(s,'%% Version %d.%d.%d',3);
    Psychtoolbox.version.major=cvv(1);
    Psychtoolbox.version.minor=cvv(2);
    Psychtoolbox.version.point=cvv(3);

    if any(strcmp(PsychtoolboxRoot, {'/usr/share/octave/site/m/psychtoolbox-3/', '/usr/share/matlab/site/m/psychtoolbox-3/'}))
        % It is a Debian version of the package
        Psychtoolbox.version.flavor = 'Debian package';
        [status, result] = system('zcat /usr/share/doc/psychtoolbox-3-common/changelog.Debian.gz| head -1 | sed -e "s/).*/)/g"');
        if status == 0
            Psychtoolbox.version.revstring = result(1:end-1);
        else
            Psychtoolbox.version.revstring = 'WARNING: failed to obtain Debian revision';
        end

        if IsOctave
            infourl = 'http://neuro.debian.net/pkgs/octave-psychtoolbox-3.html';
        else
            infourl = 'http://neuro.debian.net/pkgs/matlab-psychtoolbox-3.html';
        end
        
        % Build final version string:
        Psychtoolbox.version.string = sprintf('%d.%d.%d - Flavor: %s - %s\nFor more info visit:\n%s', Psychtoolbox.version.major, Psychtoolbox.version.minor, Psychtoolbox.version.point, ...
            Psychtoolbox.version.flavor, Psychtoolbox.version.revstring, infourl);
        
        % Retrieve the date of the Debian release:
        Psychtoolbox.date = sscanf(result, 'psychtoolbox-3 (%*d.%*d.%*d.%d.%*s');
    else
        % Additional parser code for SVN information. This is slooow!
        svncmdpath = GetSubversionPath;

        % Find revision string for Psychtoolbox that defines the SVN revision
        % to which this working copy corresponds:
        if ~IsWin
            [status , result] = system([svncmdpath 'svnversion -c ' PsychtoolboxRoot]);
        else
            [status , result] = dos([svncmdpath 'svnversion -c ' PsychtoolboxRoot]);
        end

        if status==0
            % Parse output of svnversion: Find revision number of the working copy.
            colpos=strfind(result, ':');
            if isempty(colpos)
                Psychtoolbox.version.revision=sscanf(result, '%d',1);
            else
                cvv = sscanf(result, '%d:%d',2);
                Psychtoolbox.version.revision=cvv(2);
            end
            if isempty(strfind(result, 'M'))
                Psychtoolbox.version.revstring = sprintf('Corresponds to SVN Revision %d', Psychtoolbox.version.revision);
            else
                Psychtoolbox.version.revstring = sprintf('Corresponds to SVN Revision %d but is *locally modified* !', Psychtoolbox.version.revision);
            end

            % Ok, now find the flavor and such...
            if ~IsWin
                [status , result] = system([svncmdpath 'svn info --xml ' PsychtoolboxRoot]); %#ok<*ASGLU>
            else
                [status , result] = dos([svncmdpath 'svn info --xml ' PsychtoolboxRoot]);
            end

            % First test for end-user branch:
            marker = '/github.com/Psychtoolbox-3/Psychtoolbox-3/branches/';
            startdel = strfind(result, marker) + length(marker);

            if isempty(startdel)
                % Nope: Search for developer branch aka 'trunk' aka 'master':
                marker = '/github.com/Psychtoolbox-3/Psychtoolbox-3/';
                startdel = strfind(result, marker) + length(marker);
            end

            if isempty(startdel)
                % Nope: Retry with a different query for older svn clients:
                if ~IsWin
                    [status , result] = system([svncmdpath 'svn info ' PsychtoolboxRoot]);
                else
                    [status , result] = dos([svncmdpath 'svn info ' PsychtoolboxRoot]);
                end

                % Retry first test for end-user branch:
                marker = '/github.com/Psychtoolbox-3/Psychtoolbox-3/branches/';
                startdel = strfind(result, marker) + length(marker);

                if isempty(startdel)
                    % Nope: Retry search for developer branch aka 'trunk' aka 'master':
                    marker = '/github.com/Psychtoolbox-3/Psychtoolbox-3/';
                    startdel = strfind(result, marker) + length(marker);
                end
            end

            findel = min(strfind(result(startdel:length(result)), '/Psychtoolbox')) + startdel - 2;
            Psychtoolbox.version.flavor = result(startdel:findel);

            % Retrieve the date of last commit:
            startdel = strfind(result, '<date>') + length('<date>');
            findel = strfind(result, 'T') - 1;
            Psychtoolbox.date = result(startdel:findel);

            % Build final version string:
            Psychtoolbox.version.string = sprintf('%d.%d.%d - Flavor: %s - %s\nFor more info visit:\n%s', Psychtoolbox.version.major, Psychtoolbox.version.minor, Psychtoolbox.version.point, ...
                Psychtoolbox.version.flavor, Psychtoolbox.version.revstring, Psychtoolbox.version.websvn);
        else
            % Fallback path if svn commands fail for some reason. Output as much as we can.
            fprintf('PsychtoolboxVersion: WARNING - Could not query additional version information from SVN -- svn tools not properly installed?!?\n');
            Psychtoolbox.version.string=sprintf('%d.%d.%d', Psychtoolbox.version.major, Psychtoolbox.version.minor, Psychtoolbox.version.point);
            ss=s(n:end);
            Psychtoolbox.date=ss(min(find(ss-' ')):end); %#ok<MXFND>
        end
    end
    
    if isempty(Psychtoolbox.version.websvn)
        % Build final GIT URL: This is the location where one can find detailed info about this working copy:
        % MK TODO: This doesn't have much value as long as we can't
        % correlate SVN revisions to GIT SHA-1's or something similar.
        % Just put in our main GitHub repo URL:
        Psychtoolbox.version.websvn = 'https://github.com/Psychtoolbox-3/Psychtoolbox-3';
    end
end

versionString=Psychtoolbox.version.string;
versionStructure=Psychtoolbox.version;
