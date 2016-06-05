function PsychtoolboxPostInstallRoutine(isUpdate, flavor)
% PsychtoolboxPostInstallRoutine(isUpdate [, flavor])
%
% Psychtoolbox post installation routine. You should not call this
% function directly! This routine is called by DownloadPsychtoolbox,
% or UpdatePsychtoolbox after a successfull download/update of
% Psychtoolbox. The routine performs tasks that are common to
% downloads and updates, so they can share their code/implementation.
%
% As PsychtoolboxPostInstallRoutine itself is downloaded or updated,
% it can contain code specific to each Psychtoolbox revision/release
% to perform special setup procedures for new features, to announce
% important info to the user, whatever...
%
% Currently the routine performs the following tasks:
%
% 1. Clean up the Matlab/Octave path to Psychtoolbox: Remove unneeded .svn subfolders.
% 2. Contact the Psychtoolbox server to perform online registration of this
%    working copy of Psychtoolbox.
% 3. Add the PsychJava subfolder to the static Matlab class-path if neccessary.
%    This enables the Java-based GetChar support on Matlab.
% 4. Add the PsychStartup.m routine to Matlab's startup.m file on Windows.
% 5. Perform post-installation checks and basic troubleshooting.

%
% History:
% 23/06/2006 Written (MK).
% 17/09/2006 Made working on Matlab-5 and Octave. Made more robust. (MK)
% 22/09/2006 Replace system copy commands by Matlabs copyfile() - More
%            robust (MK).
% 14/10/2006 Update web page pointers at end, just to point at new wiki (DHB).
% 28/10/2006 Accept 'current' as synonym for 'beta'. (DHB)
% 23/05/2007 Add Matlab R2007 vs. earlier detection to Windows version (MK).
% 16/04/2008 Write/Read PTB flavor to/from users PsychtoolboxConfigDir as well for higher robustness (MK).
% 15/04/2009 Add warning about unsupported OS/X systems older than Tiger (MK).
% 15/06/2009 Add support for postinstall for Octave-3.2.x, remove Octave-2 support (MK).
% 03/01/2010 Extend Matlab R2007 vs. earlier detection for Windows up to year 2014 (MK).
% 03/04/2010 Move PsychtoolboxRegistration to the end, after all real
%            installation has been finished. Add additional error-checking
%            and troubleshooting for peculiarities of Octave on Windows. (MK)
% 05/27/2010 Update instructions for downloading the vcredist_x86.exe
%            security update for users of MS-Windows. (MK)
% 12/27/2010 Add check for unsupported Matlab versions prior to V6.5. (MK)
% 10/31/2011 Add call to SwitchToNewPsychtoolboxHoster for switch to
%            GoogleCode, if needed. (MK)
% 01/06/2012 Add support for calling PsychLinuxConfiguration on Linux. (MK)
% 04/30/2012 Add support for 64-Bit OSX. (MK)
% 06/13/2012 Removed call to SwitchToNewPsychtoolboxHoster, no longer needed (DN)
% 07/10/2012 Use textscan() on R2012a+ and verLessThan() to detect R2007a+ (MK)
% 09/05/2012 Update support for 64-Bit Octave and versions > 3.2. (MK)
% 09/13/2012 Add startup.m setup for 64-Bit Matlab + 64-Bit Windows. (MK)
% 09/14/2012 Cancel support for Octave on MS-Windows. (MK)
% 09/14/2012 Cancel support for 32-Bit Octave on OSX. (MK)
% 11/11/2012 More cleanup. E.g., don't warn about Octave > 3.2 anymore. (MK)
% 04/16/2013 Use javaclasspath.txt instead of classpath.txt on R2013a and later. (MK)
% 05/13/2013 Factor out Java classpath setup to call to PsychJaveTrouble(1). (MK)
% 07/02/2013 Drop support for 32-Bit Matlab on OSX, and thereby for 32-Bit OSX. (MK)
% 07/02/2013 Drop support for OSX versions older than 10.6 "Snow Leopard". (MK)
% 09/12/2013 Setup PsychStartup.m in startup.m for 32-Bit Windows as well. (MK)
% 05/18/2014 No support for 32-Bit Matlab on Linux and Windows anymore for 3.0.12. (MK)
% 09/23/2014 No support for OSX 10.7 and earlier anymore. (MK)
% 10/05/2014 Add some request for donations at the end. (MK)
% 10/17/2015 Also add call to PsychStartup() to Octave startup for MS-Windows. (MK)
% 01/27/2016 Use Octave3 folder for mex files for both Octave-3 and Octave-4. (MK)
% 03/15/2016 Need liboctave-dev package for symlinks liboctinterp.so -> Octave specific liboctinterp.x.so (MK)
% 04/01/2016 64-Bit Octave-4 support for MS-Windows established. (MK)
% 05/13/2016 Add new rpath fixup for Octave on OSX - Copy runtime libs to search dirs. (MK)
% 06/01/2016 32-Bit Octave-4 support for MS-Windows removed. (MK)
%            64-Bit Octave-3 support for OSX removed. (MK)
%            OSX 10.8 and 10.9 support removed. (MK)

fprintf('\n\nRunning post-install routine...\n\n');

if nargin < 1
   error('PsychtoolboxPostInstallRoutine: Required argument isUpdate missing!');
end

if nargin < 2
    % No flavor provided: Default to 'unknown', but try to determine it from the
    % flavor file if this is an update.
    flavor = 'unknown';
    try
        if isUpdate>0
            % This is an update of an existing working copy. Check if flavor-file
            % is available:
            flavorfile = [PsychtoolboxRoot 'ptbflavorinfo.txt'];
            if exist(flavorfile, 'file')
                fd=fopen(flavorfile);
                if fd > -1
                    flavor = fscanf(fd, '%s');
                    fclose(fd);
                end
            end

            % Still unknown?
            if strcmp(flavor, 'unknown')
                % Yep: Retry in users PsychtoolboxConfigDir:
                flavorfile = [PsychtoolboxConfigDir 'ptbflavorinfo.txt'];
                if exist(flavorfile, 'file')
                    fd=fopen(flavorfile);
                    if fd > -1
                        flavor = fscanf(fd, '%s');
                        fclose(fd);
                    end
                end
            end
        end
    catch %#ok<*CTCH>
        fprintf('Info: Failed to determine flavor of this Psychtoolbox. Not a big deal...\n');
    end
else
    % Handle 'current' as synonym for 'beta', and 'unsupported' as synonym
    % for former 'stable'.
    switch (flavor)
        case 'current'
            flavor = 'beta';
        case 'unsupported'
            flavor = 'stable';
    end

    % Flavor provided: Write it into the flavor file for use by later update calls:
    try
        flavorfile = [PsychtoolboxRoot 'ptbflavorinfo.txt'];
        fd=fopen(flavorfile, 'wt');
        if fd > -1
            fprintf(fd, '%s\n', flavor);
            fclose(fd);
        end
    catch
        fprintf('Info: Failed to store flavor of this Psychtoolbox to file. Not a big deal...\n');
        % Retry with users PsychtoolboxConfigDir:
        try
            flavorfile = [PsychtoolboxConfigDir 'ptbflavorinfo.txt'];
            fd=fopen(flavorfile, 'wt');
            if fd > -1
                fprintf(fd, '%s\n', flavor);
                fclose(fd);
            end
        catch
            fprintf('Info: Failed to store flavor of this Psychtoolbox to file a 2nd time. Not a big deal...\n');
        end
    end
end

% Get rid of any remaining .svn folders in the path.
try
    path(RemoveSVNPaths);
    if exist('savepath') %#ok<EXIST>
        savepath;
    else
        path2rc;
    end
catch
    fprintf('Info: Failed to remove .svn subfolders from path. Not a big deal...\n');
end

% 32-Bit Octave or 32-Bit Matlab on OSX? This is unsupported as of Version 3.0.11.
if (IsOSX || IsWin) && ~Is64Bit
    fprintf('Psychtoolbox 3.0.13 and later versions do no longer work with 32-Bit versions of Octave or Matlab on OSX or Windows.\n');
    fprintf('You need to upgrade to a 64-Bit version of Octave or Matlab on these systems, which is fully supported.\n');
    fprintf('You can also use the alternate download function DownloadLegacyPsychtoolbox() to download\n');
    fprintf('an old legacy copy of Psychtoolbox-3.0.9, which did support 32-Bit Octave 3.2 on OSX, or use\n');
    fprintf('DownloadPsychtoolbox() with flavor ''Psychtoolbox-3.0.10'', which does support 32-Bit Matlab on OSX.\n');
    fprintf('DownloadPsychtoolbox() with flavor ''Psychtoolbox-3.0.12'', does support 32-Bit Octave-4 on Windows.\n');
    error('Tried to setup on 32-Bit Octave, which is no longer supported on OSX or Windows.');
end

if ~Is64Bit && ~IsOctave
    fprintf('Psychtoolbox 3.0.12 and later do no longer work with 32-Bit versions of Matlab.\n');
    fprintf('You need to upgrade to a supported 64-Bit version of Octave or Matlab. 32-Bit Octave is still\n');
    fprintf('supported on GNU/Linux.\n');
    fprintf('If you must use a legacy 32-Bit Matlab environment, you can call the function\n');
    fprintf('DownloadPsychtoolbox() with flavor ''Psychtoolbox-3.0.11'', which does support 32-Bit Matlab on Linux and Windows.\n');
    error('Tried to setup on 32-Bit Matlab, which is no longer supported.');
end

% Check if our own startup function is part of the startup file and add it,
% if it isn't already part of it. Currently we only need this for MS-Windows.
if IsWin
    % Is it already implanted? Then we ain't nothing to do:
    if ~IsPsychStartupImplantedInStartup
        % Nope. Does a proper file already exist?
        if IsOctave
            whereisit = [get_home_directory filesep '.octaverc'];
            if exist(whereisit, 'file') == 0
                whereisit = '';
            end
        else
          whereisit = which('startup.m');
        end

        if isempty(whereisit)
            % No: Create our own one.
            if IsOctave
                whereisit = [get_home_directory filesep '.octaverc'];
                fprintf('Creating a .octaverc file for Psychtoolbox at %s\n', whereisit);
            else
                whereisit = [PsychtoolboxRoot 'PsychInitialize' filesep 'startup.m'];
                fprintf('Creating a startup.m file for Psychtoolbox at %s\n', whereisit);
            end
        else
            fprintf('Adding PsychStartup() call to startup file for Psychtoolbox at %s\n', whereisit);
        end

        % whereisit points to the location of the existing or to be created
        % file. Open (or create) it in append mode:
        try
            fd = fopen(whereisit, 'a');
            fprintf(fd, '\n');
            fprintf(fd, '%% Call Psychtoolbox-3 specific startup function:\n');
            fprintf(fd, 'if exist(''PsychStartup''), PsychStartup; end;\n');
            fprintf(fd, '\n');
            fclose(fd);
        catch
            fprintf('WARNING: Failed to update or create startup file to add a call to PsychStartup()! Trouble ahead.\n');
        end
    end

    % Execute our startup function once manually, so it works already for
    % this session:
    PsychStartup;
end

% Check for operating system minor version on Mac OS/X.
if IsOSX
    % Query kernel version via system() call:
    [s, did]=system('uname -r');
    if s == 0
        % Parse string for kernel major number, then translate to OS
        % minor version by subtracting 4:
        minorver = sscanf(did, '%i') - 4;
    else
        % Failed to query: Assume we're good for now...
        minorver = inf;
    end

    % Is the operating system version < 10.10?
    if minorver < 10
        % Yes. This is MacOSX 10.9 or earlier, i.e., older than 10.10
        % Yosemite. PTB will not work on such an old system:
        fprintf('\n\n\n\n\n\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n');
        fprintf('Your operating system is Mac OS/X version 10.%i.\n\n', minorver);
        fprintf('This release of Psychtoolbox-3 is not compatible\n');
        fprintf('to OSX versions older than 10.10 "Yosemite".\n');
        fprintf('That means that almost all functionality will not work!\n\n');
        fprintf('You could download an older version of Psychtoolbox\n');
        fprintf('onto your system to get better results. See our Wiki for help.\n');
        fprintf('Better though, update your operating system to at least version\n');
        fprintf('10.11.5 or later, better the very latest OSX version.\n');
        fprintf('\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n\n');
        fprintf('Press any key on keyboard to try to continue with setup, although\n');
        fprintf('this will likely fail soon and leave you with a dysfunctional toolbox.\n\n');
        pause;
    end
end

% Matlab specific setup:
if ~IsOctave
    % Check if this is Matlab of version prior to V 7.4:
    v = ver('matlab');
    if ~isempty(v)
        v = v(1).Version; v = sscanf(v, '%i.%i.%i');
        if (v(1) < 7) || ((v(1) == 7) && (v(2) < 4))
            % Matlab version < 7.4 detected. This is no longer
            % supported.
            fprintf('\n\nYou are using a Matlab version older than Version 7.4.\n');
            fprintf('The current "beta" flavor is no longer compatible with your version of Matlab.\n');
            fprintf('Current "beta" only works on Matlab Version 7.4 (R2007a) or later.\n\n');
            fprintf('I will try to finish setup, but most functions will not work for you.\n');
            fprintf('Please run the legacy DownloadLegacyPsychtoolbox() downloader to download an outdated,\n');
            fprintf('but functional older version of Psychtoolbox (e.g., V3.0.9) for your Matlab setup or to\n');
            fprintf('receive further instructions.\n');
            fprintf('\n\nPress any key to continue after you have read and understood above message completely.\n\n');
            pause;
        end
    end
end

% Special case handling for Octave:
if IsOctave
    % GNU/Octave. Need to prepend the proper folder with
    % the pseudo-MEX files to path:
    rc = 0; %#ok<NASGU>
    rdir = '';

    try
        % Remove binary MEX folders from path:
        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFiles']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFiles64'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFiles64']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFilesARM'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3LinuxFilesARM']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles64'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles64']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4OSXFiles64'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4OSXFiles64']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3WindowsFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3WindowsFiles']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4WindowsFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4WindowsFiles']);
        end

        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4WindowsFiles64'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave4WindowsFiles64']);
        end

        % Encode prefix and Octave major version of proper folder:
        octavev = sscanf(version, '%i.%i');
        octavemajorv = octavev(1);
        octaveminorv = octavev(2);

        fprintf('Octave major version %i detected. Will prepend the following folder to your Octave path:\n', octavemajorv);

        % Octave-3 and Octave-4 can share the same mex files in the Octave-3
        % folder on Linux:
        if ismember(octavemajorv, [3,4]) && IsLinux
            rdir = [PsychtoolboxRoot 'PsychBasic' filesep 'Octave3'];
        else
            rdir = [PsychtoolboxRoot 'PsychBasic' filesep 'Octave' num2str(octavemajorv)];
        end

        % Add proper OS dependent postfix:
        if IsLinux
            rdir = [rdir 'LinuxFiles'];
        end

        if IsOSX
            rdir = [rdir 'OSXFiles'];
        end

        if IsWin
            rdir = [rdir 'WindowsFiles'];
        end

        if IsARM
            % ARM processor architecture:
            rdir = [rdir 'ARM'];
        end

        if Is64Bit
            % 64 bit Octave. Select 64 bit mex file folder:
            rdir = [rdir '64'];
        end

        fprintf(' %s ...\n', rdir);
        addpath(rdir);

        rc = savepath;
    catch
        rc = 2;
    end

    if rc > 0
        fprintf('=====================================================================\n');
        fprintf('ERROR: Failed to prepend folder %s to Octave path!\n', rdir);
        fprintf('ERROR: This will likely cause complete failure of PTB to work.\n');
        fprintf('ERROR: Please fix the problem (maybe insufficient permissions?)\n');
        fprintf('ERROR: If everything else fails, add this folder manually to the\n');
        fprintf('ERROR: top of your Octave path.\n');
        fprintf('ERROR: Trying to continue but will likely fail soon.\n');
        fprintf('=====================================================================\n\n');
    end

    if (~IsLinux && (octavemajorv < 4)) || (octavemajorv < 3) || (octavemajorv == 3 && octaveminorv < 2)
        fprintf('\n\n=================================================================================\n');
        fprintf('WARNING: Your version %s of Octave is obsolete. We strongly recommend\n', version);
        if IsLinux
            % On Linux everything >= 3.2 is fine:
            fprintf('WARNING: using the latest stable version of the Octave 3.2.x series or later for use with Psychtoolbox.\n');
        else
            % On Windows/OSX we only care about >= 4.0 atm:
            fprintf('WARNING: using the latest stable version of the Octave 4.0.x series for use with Psychtoolbox.\n');
        end
        fprintf('WARNING: Stuff may not work at all or only suboptimal with other versions and we\n');
        fprintf('WARNING: don''t provide any support for such old versions.\n');
        fprintf('\nPress any key to continue with setup.\n');
        fprintf('=================================================================================\n\n');
        pause;
    end

    if IsOSX
        % Need to copy the Octave runtime libraries somewhere our mex files can find them. The only low-maintenance
        % way of dealing with this mess of custom library pathes per octave version, revision and packaging format.
        % Preferred location is the folder with our mex files - found by rpath = @loader_path
        if ~copyfile([octave_config_info.octlibdir filesep 'liboctinterp.3.dylib'], [rdir filesep], 'f') || ...
           ~copyfile([octave_config_info.octlibdir filesep 'liboctave.3.dylib'], [rdir filesep], 'f')
            % Copy into our mex files folder failed. A second location where the linker will search is the
            % $HOME/lib directory of the current user, so try that as target location:
            tdir = PsychHomeDir('lib');
            fprintf('\n\nFailed to copy Octave runtime libraries to mex file folder [%s].\nRetrying in users private lib dir: %s ...\n', rdir, tdir);
            if ~copyfile([octave_config_info.octlibdir filesep 'liboctinterp.3.dylib'], tdir, 'f') || ...
               ~copyfile([octave_config_info.octlibdir filesep 'liboctave.3.dylib'], tdir, 'f')
                fprintf('\nFailed to copy runtime libs to [%s] as well :(.\n', tdir);
                fprintf('Our mex files will likely not work this way. Maybe the directories lack file write permissions?\n');
                fprintf('\n\n\nA last workaround would be to restart octave from a terminal via this line:\n\nexport DYLD_LIBRARY_PATH=%s ; octave\n\n\n', octave_config_info.octlibdir);
            end
        end
    end

    try
        % Rehash the Octave toolbox cache:
        path(path);
        rehash;
        clear WaitSecs;
    catch
        fprintf('WARNING: Rehashing the Octave toolbox cache failed. I may fail and recommend\n');
        fprintf('WARNING: quitting and restarting Octave, then retry.\n');
    end

    try
        % Try if WaitSecs MEX file works...
        WaitSecs(0.1);
    catch
        % Failed! Either screwed setup of path or missing runtime
        % libraries.
        fprintf('ERROR: WaitSecs-MEX does not work, most likely other MEX files will not work either.\n');
        if ismember(octavemajorv, [3,4]) && IsLinux
            fprintf('ERROR: Make sure to have the ''liboctave-dev'' package installed, otherwise symlinks\n');
            fprintf('ERROR: from liboctinterp.so to the liboctinterp library of your Octave installation\n');
            fprintf('ERROR: might by missing, causing our mex files to fail to load with linker errors.\n');
        end
        fprintf('ERROR: One reason might be that your version %s of Octave is incompatible. We recommend\n', version);
        fprintf('ERROR: use of the latest stable version of Octave-4 as announced on the www.octave.org website.\n');
        fprintf('ERROR: Another conceivable reason would be missing or incompatible required system libraries on your system.\n\n');
        fprintf('ERROR: After fixing the problem, restart this installation/update routine.\n\n');
        fprintf('\n\nInstallation aborted. Fix the reported problem and retry.\n\n');
        return;
    end

    % End of special Octave setup.
end

% Special case handling for different Matlab releases on MS-Windoze:
if IsWin && ~IsOctave
    rc = 0; %#ok<NASGU>

    if strfind(cd,'system32')
        % the below code fails if the current directory is system32 (e.g.
        % C:\Windows\system32), as it contains dlls like version.dll, which
        % get called instead of the built-in functions....
        cd(PsychtoolboxRoot);
    end

    try
        % Remove DLL folders from path:
        rmpath([PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);

        % Is this a Release2007a (Version 7.4.0) or later Matlab?
        if ~exist('verLessThan') || verLessThan('matlab', '7.4.0') %#ok<EXIST>
            % This is a pre-R2007a Matlab: No longer supported by V 3.0.10+
            fprintf('Matlab release prior to R2007a detected. This version is no longer\n');
            fprintf('supported by Psychtoolbox 3.0.10 and later. Aborted.');
            fprintf('\n\nInstallation aborted. Fix the reported problem and retry.\n\n');
            return;
        else
            % This is a R2007a or post R2007a Matlab:
            % Add PsychBasic/MatlabWindowsFilesR2007a/ subfolder to Matlab
            % path:
            rdir = [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\'];
            fprintf('Matlab release 2007a or later detected. Will prepend the following\n');
            fprintf('folder to your Matlab path: %s ...\n', rdir);
            addpath(rdir);
        end

        rc = savepath;
    catch
        rc = 2;
    end

    if rc > 0
        fprintf('=====================================================================\n');
        fprintf('ERROR: Failed to prepend folder %s to Matlab path!\n', rdir);
        fprintf('ERROR: This will likely cause complete failure of PTB to work.\n');
        fprintf('ERROR: Please fix the problem (maybe insufficient permissions?)\n');
        fprintf('ERROR: If everything else fails, add this folder manually to the\n');
        fprintf('ERROR: top of your Matlab path.\n');
        fprintf('ERROR: Trying to continue but will likely fail soon.\n');
        fprintf('=====================================================================\n\n');
    end

    try
        % Rehash the Matlab toolbox cache:
        path(path);
        rehash('pathreset');
        rehash('toolboxreset');
        clear WaitSecs;
    catch
        fprintf('WARNING: Rehashing the Matlab toolbox cache failed. I may fail and recommend\n');
        fprintf('WARNING: quitting and restarting Matlab, then retry.\n');
    end

    try
        % Try if WaitSecs MEX file works...
        WaitSecs('YieldSecs', 0.1)
    catch
        % Failed! Either screwed setup of path or missing VC++ 2010 runtime libraries.
        fprintf('ERROR: Most likely cause: The Microsoft Visual C++ 2010 runtime libraries\n');
        fprintf('ERROR: are missing on your system.\n\n');
        if IsWin(1)
            % Need 64-Bit runtime:
            fprintf('ERROR: Execute the installer file vcredist_x64.exe, which is located in your Psychtoolbox/PsychContributed/ folder.\n');
        else
            % Need 32-Bit runtime:
            fprintf('ERROR: Execute the installer file vcredist_x86.exe, which is located in your Psychtoolbox/PsychContributed/ folder.\n');
        end
        fprintf('ERROR: You must execute that installer as an administrator user. Exit Matlab before the installation, then restart it.\n');
        fprintf('ERROR: After fixing the problem, restart this installation/update routine.\n\n');
        fprintf('ERROR: You can also just do a: cd(PsychtoolboxRoot); SetupPsychtoolbox; PsychtoolboxRegistration(%i, ''%s'');\n\n', isUpdate, flavor);
        fprintf('ERROR: This will avoid a full download of Psychtoolbox over the internet and just finish the setup.\n');
        
        fprintf('\n\nInstallation aborted. Fix the reported problem and retry.\n\n');
        return;
    end
end

% If we're using Matlab then add the PsychJava stuff to the static
% Java classpath.
if ~IsOctave
    % Try to setup Matlab static Java class path:
    PsychJavaTrouble(1);
end

% Check if Screen is functional:
try
    % Linux specific instructions:
    if IsLinux
        fprintf('\n\n');
        fprintf('The Psychtoolbox on GNU/Linux needs the following 3rd party libraries\n');
        fprintf('in order to function correctly. If you get "Invalid MEX file errors",\n');
        fprintf('or similar fatal error messages, check if these are installed on your\n');
        fprintf('system and if they are missing, install them via your system specific\n');
        fprintf('software management tools:\n');
        fprintf('\n');
        fprintf('For Screen() and OpenGL support:\n\n');
        fprintf('* The OpenGL utility toolkit GLUT: glut, glut-3 or freeglut are typical provider packages in most Linux distributions.\n');
        fprintf('* GStreamer multimedia framework: At least version 1.0.0 of the core runtime and the gstreamer-base plugins.\n');
        fprintf('  For optimal performance use the latest available versions.\n');
        fprintf('  A simple way to get GStreamer at least on Ubuntu Linux is to install the "rhythmbox" or\n');
        fprintf('  "totem" multimedia-players. You may need to install additional packages to play back all\n');
        fprintf('  common audio- and video file formats. See "help GStreamer".\n');
        fprintf('* libusb-1.0 USB low-level access library.\n');
        fprintf('* libdc1394 Firewire video capture library.\n');
        fprintf('* libraw1394 Firewire low-level access library.\n');
        fprintf('\n\n');
        fprintf('For PsychKinect() (See "help InstallKinect"):\n\n');
        fprintf('* libusb-1.0 USB low-level access library.\n');
        fprintf('* libfreenect-0.5: Kinect driver library version 0.5 or later.\n');
        fprintf('\n');
        fprintf('For PsychHID() support:\n\n');
        fprintf('* libusb-1.0 USB low-level access library.\n');
        fprintf('\n\n');
        fprintf('For Eyelink():\n\n');
        fprintf('* The Eyelink core libraries from the SR-Research download website.\n');
        fprintf('\n');
        fprintf('\n');
        fprintf('If you receive an installation failure soon, then please read the output of\n');
        fprintf('"help GStreamer" first and follow the installation instructions for GStreamer\n');
        fprintf('on Linux. Psychtoolbox''s Screen() command will not work without GStreamer!\n\n');

        % Additional setup instructions for embedded/mobile devices with ARM cpu required?
        % No point anymore. Not needed for RaspberryPi, and our Nexus-7
        % support is essentially dead...
        if 0 && IsARM
            fprintf('Additionally, as this is a device with ARM processor, the helper library\n');
            fprintf('libwaffle-1.so needs to be installed in a system library folder for Screen\n');
            fprintf('to work. You can find a copy of the library in the PsychContributed/ArmArch/\n');
            fprintf('subfolder of your Psychtoolbox main folder. Rename it to libwaffle-1.so.0 during\n');
            fprintf('the copy.\n');
            fprintf('Another requirement, at least as of April 2013 and Ubuntu 13.04 for the Nexus-7,\n');
            fprintf('is that you must start octave from the command line, or via some script, like this:\n');
            fprintf('LD_PRELOAD=/usr/lib/libGLESv1_CM.so octave\n');
            fprintf('This is a workaround for a small bug in octave for Nexus-7, which would cause Screen()\n');
            fprintf('to crash shortly after opening an onscreen window.\n\n');
            fprintf('\n');
            fprintf('If you run PTB on the Nexus7 make sure you do not use the "Unity" desktop since\n');
            fprintf('this might lead to unforeseen problems. Use a desktop without 3D desktop compositor instead.\n');
            fprintf('E.g. install and use XFCE (sudo apt-get install xfce4).\n');            
        end
    end

    % Check Screen:
    AssertOpenGL;

    if IsLinux
        % Setup Desktop compositor ("Compiz") to un-redirect fullscreen windows.
        % This allows to bypass desktop composition for PTB fullscreen onscreen windows,
        % so we get the deterministic timing and high performance we want for visual
        % stimulus presentation. The command is a no-op
        PsychGPUControl('FullScreenWindowDisablesCompositor', 1);
    end

    % Try to execute online registration routine: This should be fail-safe in case
    % of no network connection.
    fprintf('\n\n');
    PsychtoolboxRegistration(isUpdate, flavor);
    fprintf('\n\n\n');

    % Tell user we're successfully done:
    fprintf('\nDone with post-installation. Psychtoolbox is ready for use.\n\n\n');

catch
    fprintf('\n\n');
    fprintf('Screen() or online registration failed to work for some reason:\n\n');
    fprintf('Check the troubleshooting instructions on our Wiki (Download section \n');
    fprintf('and FAQ section, maybe also the Bugs section).\n\n');
    fprintf('Once you manage to fix the problem (simply type ''AssertOpenGL'' to verify\n');
    fprintf('that stuff works now), you do not need to run the installation routine again,\n');
    fprintf('but can start working immediately.\n\n');
    fprintf('However, we kindly ask you to execute the following command after everything works,\n');
    fprintf('so your copy gets registered by us for statistical purpose:\n\n');
    fprintf('PsychtoolboxRegistration(%i, ''%s'');\n\n', isUpdate, flavor);
    fprintf('Thanks! Press RETURN or ENTER to confirm you read and understood the above message.\n');
    pause;
    fprintf('\n\n');
end

% Run Linux post-configuration script try-catch protected:
if IsLinux
    try
        % This script modifies/extends system configuration files
        % to optimize the system for use with Psychtoolbox:
        PsychLinuxConfiguration;
    catch
    end
end

% Some goodbye, copyright and getting started blurb...
fprintf('GENERAL LICENSING CONDITIONS AND TERMS OF USE:\n');
fprintf('----------------------------------------------\n\n');
fprintf('Almost all of the material contained in the Psychtoolbox-3 distribution\n');
fprintf('is free-software and/or open-source software under a OSI (http://www.opensource.org/)\n');
fprintf('approved license. Most material is covered by the MIT license or a MIT compatible license.\n\n');

fprintf('A few internal libraries and components are covered by other free software\n');
fprintf('licenses which we understand to be compatible with the MIT license, e.g., the GNU LGPL\n');
fprintf('license, or BSD and Apache-2.0 licenses, or they are in the public domain.\n\n');
fprintf('3rd-party components which are freely useable and redistributable for non-commercial\n');
fprintf('research use, due to the authors permissions, but are not neccessarily free / open-source\n')
fprintf('software, can be found in the "PsychContributed" subfolder of the Psychtoolbox distribution,\n');
fprintf('accompanied by their respective licenses.\n\n');

fprintf('A few components are licensed under the GNU GPL v2 license with a special linking\n');
fprintf('exception for use with Mathworks proprietary Matlab application. A very few components,\n');
fprintf('e.g., the PsychCV() function are currently covered by the standard GPL v2 license and\n');
fprintf('cannot be used with Matlab. These are only available for use with GNU/Octave.\n\n');

fprintf('Please read the license text and copyright info in the Psychtoolbox file\n');
fprintf('License.txt carefully before you use or redistribute Psychtoolbox-3.\n');
fprintf('Use of Psychtoolbox-3 components implies that you have read, understood and accepted\n');
fprintf('the licensing conditions.\n\n');
fprintf('However, in a nutshell, if you just use Psychtoolbox for your research, our licenses\n');
fprintf('don''t restrict you in any practically relevant way. Commercial users, developers or\n');
fprintf('redistributors should make sure they understood the licenses for the components they use.\n');
fprintf('If in doubt, contact one of the Psychtoolbox developers, or the original authors of the\n');
fprintf('components you want to use, modify, merge or redistribute with other software.\n\n');
fprintf('Your standard Psychtoolbox distribution comes without the source code for\n');
fprintf('the binary plugins (the MEX files). If you want to access the corresponding\n');
fprintf('source code, please type "help UseTheSource" for specific download instructions.\n\n');
fprintf('BEGINNERS READ THIS:\n');
fprintf('--------------------\n\n');
fprintf('If you are new to the Psychtoolbox, you might try this: \nhelp Psychtoolbox\n\n');
fprintf('Psychtoolbox website:\n');
fprintf('web http://www.psychtoolbox.org -browser\n');
fprintf('\n');
fprintf('Please make sure that you have a look at the PDF file Psychtoolbox3-Slides.pdf\n');
fprintf('in the Psychtoolbox/PsychDocumentation subfolder for an overview of differences\n');
fprintf('between Psychtoolbox-2 and Psychtoolbox-3 and proper use of basic features. That\n');
fprintf('folder contains various additional helpful information for use of Psychtoolbox.\n\n');
fprintf('\n');
fprintf('Please also familiarize yourself with the demos contained in the PsychDemos subfolder\n');
fprintf('and its subfolders. They show best practices for many common tasks and are generally\n');
fprintf('well documented.\n');
fprintf('\n\n');
fprintf('If you find this software useful then please consider donating some money \n');
fprintf('to support its ongoing maintenance and development. See: \n');
fprintf('\n');
fprintf('http://psychtoolbox.org/donations \n');
fprintf('\n');

fprintf('\nEnjoy!\n\n');

% Clear out everything:
if IsWin
    clear all;
end

return;
