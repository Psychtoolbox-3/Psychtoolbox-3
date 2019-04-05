function PsychJavaTrouble(installStatic)
% PsychJavaTrouble -- What to do if Java-based Psychtoolbox functions fail?
%
% You probably arrived at this help text because of an error in execution
% of one of the Psychtoolbox Java based functions, e.g., ListenChar,
% CharAvail, GetChar or FlushEvents.
%
% These are common reasons for failure of Java based Psychtoolbox
% functions:
%
% 1. You run Matlab in nojvm mode, i.e., you started Matlab via 'matlab
% -nojvm'. In that case, Matlabs Java virtual machine is disabled and so
% are all Java based Psychtoolbox functions. In that case you need to
% restart Matlab with its Java VM enabled.
%
% 2. You run GNU/Octave, which doesn't support our Java based keyboard
% functions. However, then you should not ever end here, as we use an
% alternate implementation on Octave.
%
% 3. The Psychtoolbox/PsychJava/ subfolder of your working copy of
% Psychtoolbox isn't included in Matlabs static Java classpath. The
% Psychtoolbox installer/upgrader (DownloadPsychtoolbox, SetupPsychtoolbox
% or UpdatePsychtoolbox, as well as DownloadAdditionsForLinux) usually
% tries to edit the 'classpath.txt' file, or since Matlab 8.1,
% 'javaclasspath.txt', of your Matlab installation in order to add the
% Psychtoolbox/PsychJava subfolder to Matlabs classpath. This procedure may
% fail due to insufficient access permissions on your system. You can
% verify this by entering 'type classpath.txt' at the Matlab prompt. The
% printed file should contain the path to the PsychJava folder. If it
% doesn't, you may want to edit the file yourself ('which classpath.txt'
% tells you the location of the file) or ask a system administrator to do
% it for you. After editing the file you need to restart Matlab. Instead of
% manual editing you can also call this function as PsychJavaTrouble(1); -
% This will try to automatically modify the classpath.txt or
% javaclasspath.txt file if your Matlab runs with sufficient permissions,
% e.g., administrator permissions.
%
% If you need a quick temporary fix for the problem, other than editing
% classpath.txt and restarting Matlab, then type 'PsychJavaTrouble' at the
% Matlab command prompt. The function will add the PsychJava folder to the
% dynamic classpath to immediately enable Java based Psychtoolbox
% functions. This fix is temporary however, it needs to be repeated after
% each restart of Matlab. Executing the command will also clear all
% variables and functions from Matlabs workspace (like 'clear all'), so
% adding it to experiment scripts may impair proper working of that
% scripts.
%
% 4. You didn't restart Matlab after the Psychtoolbox installer asked you
% to do so. -> Restart Matlab and retry.
%
% 5. The versions of GetCharJava bundled with Psychtoolbox are incompatible
% with the version of Java installed on your machine or bundled with your
% version of Matlab. If you have a Java SDK installed on your machine,
% Psychtoolbox will try to compile a matching version of GetCharJava. This
% should succeed on OS-X but is unlikely to work on Windows, because that
% system does not have a javac compiler installed by default.
%
% 6. Other reasons: Post to the Psychtoolbox forum and ask for help.
%
% Good luck!

% History:
% 20.09.2006 Written (MK)
% 13.05.2013 Add Matlab static classpath setup code as helper for routines
%            like PsychtoolboxPostInstallRoutine.m (MK)
% 07.07.2014 Help text cosmetic. (MK)
% 02.01.2019 Fix wrong argument order of strfind(), which caused redundant
%            classpath updates all these years! Thanks to user tibaurauer for
%            pointing this out. (MK)

% Only fix class path dynamically by default:
if nargin < 1 || isempty(installStatic)
    installStatic = 0;
end

% Dynamic fix requested, aka installStatic is 0 / false?
if ~installStatic
    % Yes:
    if IsOctave
        error('Running on GNU/Octave: Our Java functions are not supported or needed on that system.');
    end
    
    if ~exist('javaaddpath') %#ok<EXIST>
        error('Your version of Matlab does not support dynamic Java class path. Sorry!');
    end
    
    try
        fprintf('PsychJavaTrouble: Will now try to add the PsychJava folder to Matlabs dynamic classpath...\n');
        javapath = [PsychtoolboxRoot 'PsychJava'];
        javaaddpath(javapath);
        fprintf('PsychJavaTrouble: Added PsychJava folder to dynamic class path. Psychtoolbox Java commands should work now!\n');
        return;
    catch %#ok<CTCH>
        psychrethrow(psychlasterror);
        error('Failed to add PsychJava to classpath for unknown reason!');
    end
    
    return;
end

% If we reach this point, then we're asked to update the static classpath
% file, just as it would happen during PTB updates and installation:
% If we're using Matlab then add the PsychJava stuff to the static
% Java classpath.
if ~IsOctave
    try
        % Figure out the PsychJava path we need to add to the static Java
        % classpath.
        path_PsychJava = [PsychtoolboxRoot, 'PsychJava'];
        
        % Matlab 8.1 changes the rules about static java classpath. Lovely.
        if verLessThan('matlab', '8.1')
            % Legacy: Open up the classpath.txt file and find any PsychJava
            % entries.  If they exist, remove them, and put the current one
            % in the file.  This only allows on PsychJava to be on the
            % path.
            classpathFile = which('classpath.txt');
        else
            % Matlab version 8.1 (R2013a) or later. classpath.txt can't be
            % used anymore. Now they want us to store static classpath
            % definitions in a file called javaclasspath.txt inside the
            % Matlab preference folder:
            
            % Try to find the file, if it already exists, e.g., inside the
            % Matlab startup folder:
            classpathFile = which('javaclasspath.txt');
            
            % Found it?
            if isempty(classpathFile)
                % Nope. So we try the preference folder.
                % Retrieve path to preference folder. Create the folder if it
                % doesn't already exist:
                prefFolder = prefdir(1);
                classpathFile = [prefFolder filesep 'javaclasspath.txt'];
                if ~exist(classpathFile, 'file')
                    fid = fopen(classpathFile, 'w');
                    fclose(fid);
                end
            end
        end
        
        % Define name of backup file:
        bakclasspathFile = [classpathFile '.bak'];
        
        if ~verLessThan('matlab', '7.14')
            % New style method: (textread() is deprecated as of at least R2012a)
            fid = fopen(classpathFile);
            fileContentsWrapped = textscan(fid, '%s', 'delimiter', '\n');
            fclose(fid);
            fileContents = fileContentsWrapped{1};
        else
            fileContents = textread(classpathFile, '%s', 'delimiter', '\n'); %#ok<REMFF1>
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
            elseif ~isempty(strfind(fileContents{i}, 'PsychJava')) && ~pathInserted
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
                error(['Could not make a backup copy of Matlab''s JAVA path definition file. The system reports: ', w]);
            end
            madeBackup = 1; %#ok<NASGU>
            
            % Write out the new contents.
            FID = fopen(classpathFile, 'w');
            if FID == -1
                error('Could not open Matlab''s JAVA path definition file for write access.');
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
    catch %#ok<CTCH>
        lerr = psychlasterror;
        fprintf('Could not update the Matlab JAVA classpath file due to the following error:\n');
        fprintf('%s\n\n', lerr.message);
        fprintf('You likely do not have sufficient access permissions for the Matlab application\n');
        fprintf('folder or file itself to change the file %s .\n', classpathFile);
        fprintf('Please ask the system administrator to enable write-access to that file and its\n');
        fprintf('containing folder and then repeat the update procedure.\n');
        fprintf('Alternatively, ask the administrator to add the following line:\n');
        fprintf('%s\n', path_PsychJava);
        fprintf('to the file: %s\n\n', classpathFile);
        fprintf('If you skip this step, Psychtoolbox will still be mostly functional, \n');
        fprintf('but the Java-based commands ListenChar, CharAvail, GetChar and FlushEvents\n');
        fprintf('on Linux, MacOS-X and M$-Windows in Java mode will work less efficiently.\n');
        fprintf('For more info see ''help PsychJavaTrouble''.\n\n');
        fprintf('\nPress RETURN or ENTER to confirm you read and understood the above message.\n');
        pause;
        
        % Restore the old classpath file if necessary.
        if exist('madeBackup', 'var')
            [s, w] = copyfile(bakclasspathFile, classpathFile, 'f'); %#ok<NASGU,ASGLU>
        end
    end
end % if ~IsOctave
