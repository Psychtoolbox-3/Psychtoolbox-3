function callStack = AssertMex(varargin)
% Handle missing or dysfunctional MEX files on Matlab or Octave.
%
% AssertMex(targetname)
%
% AssertMex detects missing mex files.  Calling AssertMex from
% your help file, e.g. foo.m, asserts the existence of a mex file foo.mex.
% If no such mex file exists then AssertMex will exit with an error.
%
% Assert Mex is useful for detecting the error condition when a .m help
% file is mistakenly executed because the correspondig .mex file is
% missing. When foo.mex is missing, MATLAB silently executes the help file
% foo.m instead.  Calling AssertMex within foo.m detects and reports that
% error. E.g., in Screen.m you'll find AssertMex('Screen.m') to handle
% missing or dysfunctional Screen mex files properly.
%  _________________________________________________________________________
%
% See also: AssertOpenGL, COMPUTER, IsOSX, IsWin

% HISTORY:
%
% mm/dd/yy
%
% 12/14/03  awi     Wrote it
% 07/02/04  awi     Improved documentation: Added See also, explained
%                   behavior when no arguments are supplied.
%                   Changed name to AssertMex from PsychAssertMex
% 10/03/04  awi     Modified to accept 'WINDOWS', 'OS9', 'OSX'
% 01/08/05  fwc     fixed little bug, see comments
% 01/11/05  awi     Merged fwc fix into psychtoolbox.org master copy.
% 01/29/05  dgp     Cosmetic.
% 03/07/05  fwc     Fixed bugs due to which failed Assertmex to accept 'OSX', 'OS9' as valid input.
%					In line 80/81 badNames was often {''} which is not empty for isempty().
%					Also changed || into | in line 90
% 03/08/05  dgp     Add support for 'WIN'.
% 10/06/05	awi     Fixed bug: Changed "okSupNameMatches"  to match addition of  'WIN' to
%								to "okSupNames"
% 3/3/06    awi     Rewrote help for improved clarity.
% 6/13/09    mk     Update to handle execution failures and missing files
%                   on Octave-3.2 et al. as well.
% 5/28/12    mk     Update: 64-Bit OSX and Linux supported, but OSX PowerPC
%                   and Matlab < V7.4 (aka R2007a) no longer supported.

persistent okNames mexExtensions;

% Different processing on Octave build.
if IsOctave
    myName = char(varargin{1});
    if isempty(strfind(myName, '.m'))
        return;
    end
    
    octFilename = [ myName(1:end-1) 'mex'];
    fprintf('\nIn place of the expected Octave .mex binary plugin file this placeholder file was executed:\n\n');
    fprintf(['  ' myName '\n\n']);
    fprintf('This MEX file seems to be missing or inaccessible on your Octave path or it is dysfunctional:\n\n')
    fprintf(['  ' octFilename '\n\n']);
    
    fpath = which(myName(1:end-2));
    if isempty(fpath)
        fprintf('Hmm. I cannot find the file on your Octave path?!?\n\n');
        fprintf('One reason could be that your Octave path is wrong or not up to date\n');
        fprintf('for the current Psychtoolbox. You may want to run SetupPsychtoolbox to \n');
        fprintf('fix possible path problems.\n');
        fprintf('Make sure that the path is readable by you as well...\n');
    else
        % Some file with such a function exists on the Octave path:
        [dum1, dum2, fpathext] = fileparts(fpath); %#ok<ASGLU>
        % Proper extension for the .mex file?
        if isempty(strfind(fpathext, '.mex'))
            % Nope, wrong file bound:
            if Is64Bit
                oext = ['64' filesep];
            else
                oext = filesep;
            end
            
            if IsLinux
                fprintf('The following directory should be the *first one* on your Octave path:\n %s \n\n', [PsychtoolboxRoot 'PsychBasic/Octave3LinuxFiles' oext]);
            end
            if IsOSX
                fprintf('The following directory should be the *first one* on your Octave path:\n %s \n\n', [PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles' oext]);
            end
            if IsWindows
                fprintf('The following directory should be the *first one* on your Octave path:\n %s \n\n', [PsychtoolboxRoot 'PsychBasic\Octave3WindowsFiles' oext]);
            end
            
            fprintf('\n\nIf you have just installed Psychtoolbox, it worked properly and suddenly stopped working\n');
            fprintf('after you have restarted Octave, then you may encounter a bug present in Octave 3.2.0.\n');
            fprintf('Try the following at the Octave prompt: First type "savepath" + Enter to save the current Octave path again.\n');
            fprintf('Then exit Octave and restart it. After that try again - It may work then, if the problem was caused by aforementioned bug.\n\n');
        else
            % Correct file with correct extension, still load failure:
            % Check for supported Octave version:
            curversion = sscanf(version, '%i.%i.%i');
            if curversion(1) < 3 || curversion(2) < 2
                fprintf('Your version of Octave (%s) is incompatible with Psychtoolbox: We support Octave 3.2.0 or later.\n', version);
                error('Tried to run Psychtoolbox on an incompatible Octave version.\n');
            end
            
            fprintf('A reason could be some missing 3rd party dynamic link shared libraries on your system.\n');
            fprintf('Our default installation also only supports 32 bit versions of operating system and Octave.\n');
            fprintf('Another reason could be some binary incompatibility. You would need to recompile Psychtoolbox from source!\n\n');
        end
    end
    error('Missing, inaccessible or dysfunctional Psychtoolbox Mex file for this system, or Octave path wrong. Read the help text above carefully!!\n');
end

% Initialize the persistent variables.
if isempty(okNames)
    okNames = {      'PCWIN',  'PCWIN64', 'SOL2', 'HPUX', 'HP700', 'ALPHA', 'IBM_RS', 'SGI', 'GLNX86', 'MAC',    'MAC2', 'MACI',    'i486-pc-linux-gnu', 'MACI64',    'GLNXA64' };
    mexExtensions = {'mexw32', 'mexw64',  '*',    '*',    '*',     '*',     '*',      '*',   'mexglx', 'mexmac', 'mex',  'mexmaci', 'mex',               'mexmaci64', 'mexa64' };
end

inputNames = [];

% Check to see if there should be a mex file for our platform.
if isempty(inputNames) || ismember(computer, inputNames)
    % Element 1 will always be AssertMex. Element 2 will be the calling
    % function unless it is invoked from the commnand line.
    callStack = dbstack;
    
    if length(callStack) > 1
        callerName = callStack(2).name;
    else
        error('AssertMex was invoked from the command line.');
    end
    
    % Generate error strings
    extensionNameIndex=find(streq(computer,okNames));
    extensionName=mexExtensions{extensionNameIndex}; %#ok<FNDSB>
    mExtensionIndices=strfind('.m',callerName);
    if ~isempty(mExtensionIndices)
        termExtension=mExtensionIndices(end);
        mexFilenameHead=callerName(1:termExtension-1);
    else
        mexFilenameHead = callerName;
    end
    mexFilename=[mexFilenameHead '.' extensionName];
    fprintf('\nIn place of the expected mex file this placeholder file was executed:\n\n');
    fprintf(['  ' callerName '\n\n']);
    fprintf('This mex file seems to be missing or inaccessible on your Matlab path or it is dysfunctional:\n\n')
    fprintf(['  ' mexFilename '\n\n']);
    
    if isempty(which(mexFilename))
        if strcmp(computer, 'MAC')
            % Mac PowerPC:
            fprintf('Matlab running on MacOSX on a PowerPC based Apple Macintosh computer is no longer supported.\n');
            fprintf('You would need to get a version 3.0.9 or earlier Psychtoolbox for this to work on your computer.\n');
            fprintf('Or upgrade your computer to run a modern version of GNU/Linux and use the Psychtoolbox for Linux on PowerPC.\n\n');
        elseif (IsWin && exist('mexext') && strcmp(mexext, 'dll')) %#ok<EXIST>
            % Windows with Matlab pre-R2007a:
            fprintf('Matlab versions older than V7.4 aka R2007a running on MS-Windows are no longer supported by this release.\n');
            fprintf('Install GNU/Octave or a more modern Matlab version instead. Or downgrade to Psychtoolbox version 3.0.9 or earlier.\n\n');
        else
            fprintf('Hmm. I cannot find the file on your Matlab path?!?\n\n');
            fprintf('One reason could be that your Matlab path is wrong or not up to date\n');
            fprintf('for the current Psychtoolbox. You may want to run SetupPsychtoolbox to \n');
            fprintf('fix possible path problems.\n\n');
        end
    else
        if IsWin
            fprintf('It is important that the folder which contains the Screen mex file is located *before*\n');
            fprintf('the PsychBasic folder on your Matlab path.\n');
            fprintf('On Matlab V7.4 (R2007a) or later versions, the folder\n');
            fprintf('%sPsychBasic\\MatlabWindowsFilesR2007a\\ must be before the folder\n%sPsychBasic\\ \n\n', PsychtoolboxRoot, PsychtoolboxRoot);
            fprintf('type path to display the current path and check for this condition.\nIf the order is wrong, ');
            fprintf('simply cd into your Psychtoolbox root folder\n%s and then run SetupPsychtoolbox again.\n', PsychtoolboxRoot);
            fprintf('That is the simplest way of fixing such path problems - Or to get more diagnostic output.\n\n');
        end
        
        fprintf('Another reason could be insufficient access permissions or \n');
        fprintf('some missing 3rd party libraries on your system.\n\n');
        
        if IsWin
            fprintf('On Microsoft Windows with supported Matlab versions (>= V7.4) it could also be that\n');
            fprintf('the required Visual C++ 2010 runtime libraries are missing on your system.\n');
            fprintf('The Psychtoolbox/PsychContributed/ subfolder contains installer files for them, which\n');
            fprintf('you can execute after quitting Matlab. The name of the file is ');
            if IsWin(1)
                fprintf('vcredist_x64.exe\n');
            else
                fprintf('vcredist_x86.exe\n');
            end
        end
    end
    error('Missing or dysfunctional Psychtoolbox Mex file for this operating system. Read the help text above carefully!!');
end
