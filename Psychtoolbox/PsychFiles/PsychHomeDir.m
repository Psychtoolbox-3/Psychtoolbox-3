function ThePath=PsychHomeDir(subDir)
% Purpose: Look for the home directory of the user.
%
% Syntax: path=PsychHomeDir([subDir])
%
%          When called without optional 'subDir' argument, the path to the
%          users home folder is returned. When 'subDir' is given, the
%          path to the subfolder 'subDir' inside the home folder
%          is returned - and the 'subDir' created inside that folder
%          if neccessary.
%
% History: 1/23/08    mpr configured it was about time to write this
%          3/7/08     mpr streamlined this
%          3/8/08     mk  A bit more of streamlining - Don't write the
%                         PsychPrefsfolder.m file anymore.
%          4/28/08    mk  Made compatible with Octave, added 'subDir'
%                         option.
%          2/06/09    mk  Derived from PsychtoolboxConfigDir().

persistent PTBHomePath %#ok<REDEF>
ThePath = [];

% Already have a cached path to config directory?
if ~isempty(PTBHomePath) %#ok<NODEF>
    if exist(PTBHomePath,'dir') %#ok<NODEF>
        % Yes - Assign it:
        ThePath=PTBHomePath;
    end
end

% No path yet? If so find it - and create configdir folder if neccessary.
if isempty(ThePath)
    if IsOSX
        % Did this instead of '~/' because the which command above and the addpath
        % commands below will expand '~/' to a full path; echoing the HOME
        % environment variable was the first way I found to get said full path so
        % that strings will match when they should
        [ErrMsg,HomeDir] = unix('echo $HOME');
        % end-1 to trim trailing carriage return
        StringStart = [HomeDir(1:(end-1)) ''];
    elseif IsLinux
        [ErrMsg,HomeDir] = unix('echo $HOME');
        % end-1 to trim trailing carriage return
        StringStart = [HomeDir(1:(end-1)) ''];
    elseif IsWindows
        [ErrMsg,StringStart] = dos('echo %UserData%');
        % end-1 to trim trailing carriage return
        %StringStart = StringStart(1:(end-1));
		StringStart = deblank(StringStart);
        if strfind(StringStart,'%UserData')
            FoundHomeDir = 0;
            [ErrMsg,HomeDir] = dos('echo %UserProfile%');
            % HomeDir = HomeDir(1:(end-1));
			HomeDir = deblank(HomeDir);
            if strfind(HomeDir,'%UserProfile')
                HomeDir = uigetdir('','Please find your home folder for me');
                if ischar(HomeDir)
                    FoundHomeDir = 1;
                else
                    warning(sprintf(['I could not find your home directory or understand your input so I am storing\n' ...
                        'preferences folder in the current working directory: %s.\n'],pwd)); %#ok<SPWRN>
                    StringStart = [pwd filesep];
                end
            else
                FoundHomeDir = 1;
            end
            if FoundHomeDir
                [DirMade,DirMessage]=mkdir(HomeDir,'Application Data'); %#ok<NASGU>
                if DirMade
                    StringStart = [HomeDir filesep 'Application Data' filesep];
                else
                    warning(sprintf('"Application Data" folder neither exists nor is createable;\nstoring preferences in home directory.')); %#ok<WNTAG,SPWRN>
                    StringStart = [HomeDir filesep];
                end
            end
        else
            StringStart = [StringStart filesep];
        end
    else
        fprintf(['I do not know your operating system, so I don''t know where I should store\n' ...
            'Preferences.  I''m putting them in the current working directory:\n      %s.\n\n'],pwd);
        StringStart = [pwd filesep];
    end

    TheDir = strtrim(StringStart);

    if exist(TheDir,'dir')
        ThePath = TheDir; %#ok<NASGU>
    else
        error(sprintf('I could not find your home directory as expected in\n\n%s\n\nWhat are the permissions on that folder?',TheDir)); %#ok<SPERR>
    end

    ThePath = [ThePath filesep];
    PTBHomePath = ThePath;
end

% Ok, Psychtoolbox root configuration folder exists and 'ThePath' is the
% fully qualified path to it.

% Did user specify a subDir inside that folder?
if exist('subDir', 'var')
    % Yes. Usercode wants path to subdirectory inside config dir. Assemble
    % path name:
    ThePath = [ThePath subDir];
    if ThePath(end) ~= filesep
        ThePath = [ThePath filesep];
    end
    
    % Create subDir on first use:
    if ~exist(ThePath, 'dir')
        [DirMade, DirMessage] = mkdir(ThePath);
        if DirMade == 0
            error(sprintf('I could not create a subfolder in your homedirectory in\n\n%s [%s]\n\nWhat are the permissions on that folder?',StringStart, DirMessage)); %#ok<SPERR>
        end
    end
end

% Return 'ThePath':
return;
