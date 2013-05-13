function ThePath = DaqtoolboxConfigDir
% Syntax: AbsolutePath = DaqtoolboxConfigDir
%
% Purpose: Look for (create if necessary) folder containing preferences for Daq
%          toolbox;
%
% History: 1/28/08  mpr configured this was needed
%          3/7/08   mpr streamlined this
%
% Function does not assume that user is using Psychophysics toolbox, but will
% subordinate Daqtoolbox if that assumption is correct.  That is, function will
% look for Daqtoolbox preferences in a folder containing Psychtoolbox
% preferences.  It will create its own folder iff it does not find one, and that
% folder will be in the Psychtoolbox preferences folder if that folder exists.

persistent DTBPrefPath

if ~isempty(DTBPrefPath)
  if exist(DTBPrefPath,'dir')
    ThePath = DTBPrefPath;
    return;
  end
end

if exist('PsychtoolboxConfigDir','file') == 2
  StringStart = [PsychtoolboxConfigDir filesep]; %#ok<NASGU>
  UsePsych=1;
else
  UsePsych=0;
  if IsOSX
    [ErrMsg,HomeDir] = unix('echo $HOME');
    % end-1 to trim trailing carriage return
    StringStart = [HomeDir(1:(end-1)) '/Library/Preferences/']; %#ok<NASGU>
  elseif IsLinux
    [ErrMsg,HomeDir] = unix('echo $HOME');    
    % end-1 to trim trailing carriage return
    StringStart = [HomeDir(1:(end-1)) '/.'];
  elseif IsWindows
    [ErrMsg,StringStart] = dos('echo %AppData%');
    % end-1 to trim trailing carriage return
    StringStart = StringStart(1:(end-1));
    if strcmp(StringStart,'%AppData%')
      FoundHomeDir = 0;
      [ErrMsg,HomeDir] = dos('echo %UserProfile%');
      HomeDir = HomeDir(1:(end-1));
      if strcmp(HomeDir,'%UserProfile%')
        HomeDir = uigetdir('','Please find your home folder for me');
        if ischar(HomeDir)
          FoundHomeDir = 1;
        else
          warning(sprintf(['I could not find your home directory or understand your input so I am storing\n' ...
                           'preferences folder in the current working directory: %s.\n'],pwd)); %#ok<WNTAG,SPWRN>
          StringStart = [pwd filesep]; %#ok<NASGU>
        end
      else
        FoundHomeDir = 1;        
      end
      if FoundHomeDir
        [DirMade,DirMessage]=mkdir(HomeDir,'Application Data'); %#ok<NASGU>
        if DirMade
          StringStart = [HomeDir filesep 'Application Data' filesep]; %#ok<NASGU>
        else
          warning(sprintf('"Application Data" folder neither exists nor is createable;\nstoring preferences in home directory.')); %#ok<WNTAG,SPWRN>
          StringStart = [HomeDir filesep]; %#ok<NASGU>
        end
      end
    else
      StringStart = [StringStart filesep];
    end
  else % if IsOSX; else
    fprintf(['I do not know your operating system, so I don''t know where I should store\n' ...
            'Preferences.  I''m putting them in the current working directory:\n      %s.\n\n'],pwd);
    StringStart = [pwd filesep];
  end % if IsOSX; else
end % if UsePsych; else
  
TheDir = [ StringStart 'Daqtoolbox'];

if exist(TheDir,'dir')
  ThePath = TheDir;
else % if exist(TheDir,'dir')
  % Could just use:
  % [DirMade, DirMessage] = mkdir(TheDir);
  % but I've had problems (in OS 9) having mkdir fail when string passed is too
  % long.  I don't know if that has been fixed in other OS's, so this is safe:
  if IsLinux && ~UsePsych
    [DirMade, DirMessage] = mkdir(StringStart(1:(end-1)),'.Daqtoolbox');
  else      
    [DirMade, DirMessage] = mkdir(StringStart,'Daqtoolbox');
  end
  if DirMade
    ThePath = [StringStart 'Daqtoolbox'];
  else 
    error(sprintf('I could not create a folder to store your preferences in\n\n%s\n\nWhat are the permissions on that folder?',StringStart)); %#ok<SPERR>
  end % if DirMade; else
end % if exist(TheDir,'dir')

DTBPrefPath = ThePath;

return;
