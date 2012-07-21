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

fprintf('\n\nRunning post-install routine...\n\n');

if nargin < 1
   error('PsychtoolboxPostInstallRoutine: Required argument isUpdate missing!');
end;

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

% Check for operating system minor version on Mac OS/X when running under
% Matlab:
if IsOSX
    if ~IsOctave & ~IsOSX(1) 
        % Running on Matlab + OS/X. Find the operating system minor version,
        % i.e., the 'y' in the x.y.z number, e.g., y=3 for 10.3.7:

        % Get 32-digit binary encoded minor version from Gestalt() MEX file:
        binminor = Gestalt('sys2');

        % Decode into decimal digit:
        minorver = 0;
        for i=1:32
            minorver = minorver + binminor(i) * 2^(32-i);
        end
    else
        % Running on Octave + OS/X: Query kernel version via system() call:
        [s, did]=system('uname -r');
        if s == 0
            % Parse string for kernel major number, then translate to OS
            % minor version by subtracting 4:
            minorver = sscanf(did, '%i') - 4;
        else
            % Failed to query: Assume we're good for now...
            minorver = inf;
        end
    end
    
    % Is the operating system minor version 'minorver' < 4?
    if minorver < 4
        % Yes. This is MacOS/X 10.3 or earlier, i.e., older than 10.4
        % Tiger. In all likelihood, this current PTB release won't work on
        % such a system anymore, because some of the binary MEX files are
        % linked against incompatible runtimes and frameworks. Output a
        % clear warning message about this, with tips on how to resolve the
        % problem:
        fprintf('\n\n\n\n\n\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n');
        fprintf('Your operating system is Mac OS/X version 10.%i.\n\n', minorver);
        fprintf('This release of Psychtoolbox-3 is likely no longer compatible\n');
        fprintf('to OS/X versions older than 10.4 "Tiger".\n\n');
        fprintf('That means that some or many crucial functions will fail.\n');
        fprintf('You may encounter errors or failures during the remainder of\n');
        fprintf('this installation procedure, or later on during use of the toolkit.\n\n');
        fprintf('You can fix this problem by switching to the last Psychtoolbox version\n');
        fprintf('that was known to (mostly) work on your old operating system:\n\n');
        fprintf('1. Type clear all on the Matlab command prompt.\n\n');
        fprintf('2. Run the downloader script DownloadPsychtoolbox again,\n');
        fprintf('   (see help DownloadPsychtoolbox for usage information),\n');
        fprintf('   but specify the optional "flavor" argument as:\n\n');
        fprintf('   ''Psychtoolbox-3.0.8-PreTiger''\n\n');
        fprintf('   Example for a standard installation in standard location:\n\n');
        fprintf('   DownloadPsychtoolbox([],[],[],''Psychtoolbox-3.0.8-PreTiger'');\n\n\n');
        fprintf('This will delete your current copy of Psychtoolbox-3 and replace it\n');
        fprintf('by the last known good version for your system.\n\n');
        fprintf('Please note that we will no longer provide support, bug fixes or enhancements\n');
        fprintf('for this old release on your old and obsolete operating system - You are on your own.\n\n');
        fprintf('We strongly recommend that you upgrade your system to a more recent OS/X version soon.\n\n');
        fprintf('Thanks for your attention and good luck!');
        fprintf('\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n\n');
        fprintf('Press any key on keyboard to continue with setup...\n');
        pause;
    end
    
    % Is the operating system minor version 'minorver' < 5 on 64-Bit OSX?
    if (minorver < 5) & IsOSX(1)
        % Yes. This is MacOS/X 10.4 or earlier, i.e., older than 10.5
        % Leopard. In all likelihood, this current PTB release won't work on
        % such a system anymore, because the binary 64-Bit MEX files are
        % linked against incompatible runtimes and frameworks. Output a
        % clear warning message about this, with tips on how to resolve the
        % problem:
        fprintf('\n\n\n\n\n\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n');
        fprintf('Your operating system is Mac OS/X version 10.%i.\n\n', minorver);
        fprintf('This 64-Bit release of Psychtoolbox-3 is not compatible\n');
        fprintf('to OS/X versions older than 10.5 "Leopard".\n\n');
        fprintf('That means that some or many crucial functions will fail.\n');
        fprintf('You may encounter errors or failures during the remainder of\n');
        fprintf('this installation procedure, or later during use of the toolkit.\n\n');
        fprintf('You can either use a version of 32-Bit Matlab to use the 32-Bit Psychtoolbox\n');
        fprintf('on your system, or a 32-Bit version of GNU/Octave. Alternatively update your\n');
        fprintf('operating system to at least version 10.5, but better 10.6 Snow Leopard or later.\n');
        fprintf('\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n\n');
        fprintf('Press any key on keyboard to continue with setup...\n');
        pause;
    end

    % Is the operating system minor version 'minorver' < 6 on 64-Bit OSX?
    if (minorver < 6) & IsOSX(1) %#ok<*AND2>
        % Yes. This is MacOS/X 10.5 or earlier, i.e., older than 10.6
        % Snow Leopard. 64-Bit PTB will only provide limited functionality:
        fprintf('\n\n\n\n\n\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n');
        fprintf('Your operating system is Mac OS/X version 10.%i.\n\n', minorver);
        fprintf('This 64-Bit release of Psychtoolbox-3 is not fully compatible\n');
        fprintf('to OS/X versions older than 10.6 "Snow Leopard".\n\n');
        fprintf('That means that some functionality will be limited or missing.\n');
        fprintf('E.g., video capture, video recording, movie playback and movie writing\n');
        fprintf('will be disabled, unless you manage to find a version of GStreamer compatible\n');
        fprintf('with OSX 10.5. High precision framebuffer display modes and some parts of\n');
        fprintf('OpenGL 3D graphics support will be disabled or limited in performance and functionality.\n');
        fprintf('General graphics performance may be lower and resource consumption may be higher.\n');
        fprintf('Windowed display mode on multi-display setups may behave slightly odd.\n\n');
        fprintf('You can either live with these restrictions, or use a version of 32-Bit Matlab or Octave\n');
        fprintf('together with the 32-Bit OSX Psychtoolbox on your system to get better results.\n');
        fprintf('Alternatively, update your operating system to at least version 10.6 "Snow Leopard" or later.\n');
        fprintf('\n\n\n==== WARNING WARNING WARNING WARNING ====\n\n\n');
        fprintf('Press any key on keyboard to continue with setup...\n');
        pause;
    end
end

% Matlab specific setup:
if ~IsOctave
    % Check if this is Matlab of version prior to V 7.4:
    v = ver('matlab');
    if ~isempty(v)
        v = v(1).Version; v = sscanf(v, '%i.%i.%i');
        if (v(1) < 7) | ((v(1) == 7) & (v(2) < 4)) %#ok<AND2,OR2>
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
    % OS/X or Linux under Octave. Need to prepend the proper folder with
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
        
        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles']);
        end
        
        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles64'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3OSXFiles64']);
        end
        
        if exist([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3WindowsFiles'], 'dir')
            rmpath([PsychtoolboxRoot 'PsychBasic' filesep 'Octave3WindowsFiles']);
        end
        
        % Encode prefix and Octave major version of proper folder:
        octavev = sscanf(version, '%i.%i');
        octavemajorv = octavev(1);
        octaveminorv = octavev(2);
        
        rdir = [PsychtoolboxRoot 'PsychBasic' filesep 'Octave' num2str(octavemajorv)];
        
        % Add proper OS dependent postfix:
        if IsLinux
            rdir = [rdir 'LinuxFiles'];
            if IsLinux(1)
                % 64 bit Octave on 64 bit Linux. Select 64 bit mex file folder:
                rdir = [rdir '64'];
            end
        end
        
        if IsOSX
            rdir = [rdir 'OSXFiles'];
            if IsOSX(1)
                % 64 bit Octave on 64 bit OSX. Select 64 bit mex file folder:
                rdir = [rdir '64'];
            end
        end
        
        if IsWin
            rdir = [rdir 'WindowsFiles'];
        end

        fprintf('Octave major version %i detected. Will prepend the following folder to your Octave path:\n', octavemajorv);
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
    
    if octavemajorv < 3 | octaveminorv < 2 %#ok<OR2>
        fprintf('\n\n=================================================================================\n');
        fprintf('WARNING: Your version %s of Octave is obsolete. We strongly recommend\n', version);
        fprintf('WARNING: using the latest stable version of the Octave 3.2.x series for use with Psychtoolbox.\n');
        fprintf('WARNING: Stuff may not work at all or only suboptimal with earlier versions and we\n');
        fprintf('WARNING: don''t provide any support for such old versions.\n');
        fprintf('\nPress any key to continue with setup.\n');
        fprintf('=================================================================================\n\n');
        pause;
    end
    
    if octavemajorv > 3 | (octavemajorv == 3 & octaveminorv > 2) %#ok<AND2,OR2>
        fprintf('\n\n=================================================================================\n');
        fprintf('INFO: Your version %s of Octave is version 3.3 or later.\n', version);
        fprintf('INFO: Psychtoolbox seems to work correctly on Octave 3.4 - 3.6, but no extensive\n');
        fprintf('INFO: systematic testing has been performed yet by the core developers.\n');
        fprintf('\nPress any key to continue with setup.\n');
        fprintf('=================================================================================\n\n');
        pause;
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
        % Try if Screen MEX file works...
        WaitSecs(0.1);
    catch
        % Failed! Either screwed setup of path or missing runtime
        % libraries.
        fprintf('ERROR: WaitSecs-MEX does not work, most likely other MEX files will not work either.\n');
        fprintf('ERROR: One reason might be that your version %s of Octave is incompatible. We recommend\n', version);        
        fprintf('ERROR: use of the latest stable version of Octave-3.2.x as announced on www.octave.org website.\n');
        fprintf('ERROR: Another conceivable reason would be missing or incompatible required system libraries on your system.\n\n');
        fprintf('ERROR: After fixing the problem, restart this installation/update routine.\n\n');
        fprintf('\n\nInstallation aborted. Fix the reported problem and retry.\n\n');
        return;
    end
    
    % End of special Octave setup.
end

% Special case handling for different Matlab releases on MS-Windoze:
if IsWin & ~IsOctave %#ok<AND2>
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
            fprintf('Execute the installer file vcredist_x64.exe, which is located in your Psychtoolbox/PsychContributed/ folder.\n');
        else
            % Need 32-Bit runtime:
            fprintf('Execute the installer file vcredist_x86.exe, which is located in your Psychtoolbox/PsychContributed/ folder.\n');
        end
        fprintf('You must execute that installer as an administrator user. Exit Matlab before the installation, then restart it.\n');
        %fprintf('Go to the following URL:\n\n');
        %fprintf('http://www.microsoft.com/downloads/details.aspx?familyid=766A6AF7-EC73-40FF-B072-9112BAB119C2&displaylang=en#filelist\n\n');
        %fprintf('ERROR: Download and install the required runtime libraries.\n\n');
        %fprintf('ERROR: Use the download button right to vcredist_x86.exe - The file with a size of 2.6 MB.\n');
        %fprintf('ERROR: Then double-click and run the downloaded vcredist_x86.exe installer to update your system.\n');
        %fprintf('ERROR: If you install the wrong runtime, it will still not work.\n\n');
        fprintf('ERROR: After fixing the problem, restart this installation/update routine.\n\n');
        fprintf('ERROR: You can also just do a: cd(PsychtoolboxRoot); SetupPsychtoolbox;\n\n');
        fprintf('ERROR: This will avoid a full download of Psychtoolbox over the internet and just finish the setup.\n');
        
        fprintf('\n\nInstallation aborted. Fix the reported problem and retry.\n\n');
        return;
    end
end

% If we're using Matlab then add the PsychJava stuff to the static
% Java classpath.
if ~IsOctave
    try
       % Figure out the PsychJava path we need to add to the static Java
        % classpath.
        path_PsychJava = [PsychtoolboxRoot, 'PsychJava'];

        % Open up the classpath.txt file and find any PsychJava entries.  If
        % they exist, remove them, and put the current one in the file.  This
        % only allows on PsychJava to be on the path.
        classpathFile = which('classpath.txt');
        bakclasspathFile = [classpathFile '.bak'];
        
        if ~verLessThan('matlab', '7.14')
            % New style method: (textread() is deprecated as of at least R2012a)
            fid = fopen(classpathFile);
            fileContentsWrapped = textscan(fid, '%s');
            fclose(fid);
            fileContents = fileContentsWrapped{1};
        else
            fileContents = textread(classpathFile, '%s'); %#ok<REMFF1>
        end
        j = 1;
        newFileContents = {};
        pathInserted = 0;
        for i = 1:length(fileContents)
            % Look for the first instance of PsychJava in the classpath and
            % replace it with the new one.  All other instances will be
            % ignored.
            if isempty(strfind(fileContents{i}, 'PsychJava'))
                newFileContents{j, 1} = fileContents{i}; %#ok<AGROW>
                j = j + 1;
            elseif ~isempty(strfind(fileContents{i}, 'PsychJava')) & ~pathInserted %#ok<AND2>
                newFileContents{j, 1} = path_PsychJava; %#ok<AGROW>
                pathInserted = 1;
                j = j + 1;
            end
        end

        % If the PsychJava path wasn't inserted, then this must be a new
        % installation, so we append it to the classpath.
        if ~pathInserted
            newFileContents{end + 1, 1} = path_PsychJava;
        end

        % Now compare to see if the new and old classpath are the same.  If
        % they are, then there's no need to do anything.
        updateClasspath = 1;
        if length(fileContents) == length(newFileContents)
            if strcmp(fileContents, newFileContents)
                updateClasspath = 0;
            end
        end

        if updateClasspath
            % Make a backup of the old classpath.
            clear madeBackup;

            [s, w] = copyfile(classpathFile, bakclasspathFile, 'f');

            if s==0
                error(['Could not make a backup copy of Matlab''s JAVA path definition file ''classpath.txt''. ' ...
                    'The system reports: ', w]);
            end
            madeBackup = 1; %#ok<NASGU>

            % Write out the new contents.
            FID = fopen(classpathFile, 'w');
            if FID == -1
                error('Could not open Matlab''s JAVA path definition file ''classpath.txt'' for write access.');
            end
            for i = 1:length(newFileContents)
                fprintf(FID, '%s\n', newFileContents{i});
            end
            fclose(FID);

            fprintf('\n\n');
            disp('*** Matlab''s Static Java classpath definition file modified. You will have to restart Matlab to enable use of the new Java components. ***');
            fprintf('\nPress RETURN or ENTER to confirm you read and understood the above message.\n');
            pause;
        end
    catch
        lerr = psychlasterror;
        fprintf('Could not update the Matlab JAVA classpath.txt file due to the following error:\n');
        fprintf('%s\n\n', lerr.message);
        fprintf('Probably you do not have sufficient access permissions for the Matlab application folder\n');
        fprintf('or the file itself to change the file %s .\n\n', classpathFile);
        fprintf('Please ask the system administrator to enable write-access to that file and its\n');        
        fprintf('containing folder and then repeat the update procedure.\n');
        fprintf('Alternatively, ask the administrator to add the following line:\n\n');
        fprintf('%s\n\n', path_PsychJava);
        fprintf('to the file: %s\n\n', classpathFile);        
        fprintf('If you skip this step, Psychtoolbox will still be mostly functional, \n');
        fprintf('but the Java-based commands ListenChar, CharAvail, GetChar and FlushEvents\n');
        fprintf('on Linux, MacOS-X and MS-Windows in Java mode will not work well - or at all.\n');
        fprintf('For more info see ''help PsychJavaTrouble''.\n\n');
        fprintf('\nPress RETURN or ENTER to confirm you read and understood the above message.\n');
        pause;

        % Restore the old classpath file if necessary.
        if exist('madeBackup', 'var')
            [s, w] = copyfile(bakclasspathFile, classpathFile, 'f'); %#ok<*ASGLU,NASGU>
        end
    end
end % if ~IsOctave

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
        fprintf('* GStreamer multimedia framework: At least version 0.10.24 of the core runtime and the gstreamer-base plugins.\n');
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
        fprintf('* libfreenect: Kinect driver library.\n');
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
    if IsOctave & IsWin %#ok<AND2>
        % Probably GStreamer runtime dll's missing:
        fprintf('Screen() or online registration failed to work under MS-Windows with GNU/Octave-3:\n\n');
        fprintf('Probably the required GStreamer multimedia framework is not yet installed on your system.\n\n');
        fprintf('Please type ''help GStreamer'' and follow the displayed installation instructions carefully.\n');
        fprintf('After this one-time setup, the Screen command should work properly.\n\n');
        fprintf('If this has been ruled out as a reason for failure, check the troubleshooting instructions on\n');
        fprintf('our Wiki (Download section and FAQ section, maybe also the Bugs section).\n\n');
    else
        fprintf('Screen() or online registration failed to work for some reason:\n\n');
        fprintf('Check the troubleshooting instructions on our Wiki (Download section \n');
        fprintf('and FAQ section, maybe also the Bugs section).\n\n');
    end
    
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

fprintf('\nEnjoy!\n\n');

% Clear out everything:
if ~IsOctave & IsWin %#ok<AND2>
    clear all;
end

return;
