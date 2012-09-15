function isRunning=IsUpdateRunning

% isRunning=IsUpdateRunning
%
% OS X: ___________________________________________________________________
%
% Report the number of processes named "update" which are running.  Only
% one update process at most should run but IsUpdateRunning will report
% higher values if indeed multiple update processes are running.
%
% Priority and Rush shut down the OS X update process when raising priority
% above 0. They use IsUpdateRunning to verify expectations about whether
% update is running prior to taking action.  
%
% OS 9: ___________________________________________________________________
%
% IsUpdateRunning does not exist in OS 9. 
%
% WINDOWS: ________________________________________________________________
% 
% IsUpdateRunning does not exist in Windows.
% 
% _________________________________________________________________________
%
% see also: GetProcessList, GetRawProcesslist, Priority, Rush, finish.

% HISTORY
% 12/06/03  awi     wrote it.
% 7/15/03   awi     wrote it.
% 3/10/09    mk     Prevent check and always return "true" on OS/X >=
%                   10.4.7 to prevent trouble with Matlab R2009a et al.

% Any need for this killupdate schnickschnack at all?
persistent killUpdateNotNeeded;
if isempty(killUpdateNotNeeded)
    % Check if this is MacOS-X 10.4.7 or later. We don't need to kill
    % the update process anymore if that is the case.
    c = Screen('Computer');
    osrelease = sscanf(c.kern.osrelease, '%i.%i.%i');

    if (osrelease(1)==8 && osrelease(2)>=7) || (osrelease(1)>=9)
        % OS-X 10.4.7 or later -> No need to kill update.
        killUpdateNotNeeded = 1;
    else
        % Pre 10.4.7 system -> Play safe and kill update.
        killUpdateNotNeeded = 0;
    end

    % Override for the scared.
    if exist('AlwaysKillUpdate', 'file')>0
        % Detected a veto file created by the user. We do kill update.
        killUpdateNotNeeded = 0;
    end
end

if killUpdateNotNeeded
    % No. We don't kill update, so we don't care about its state and simply
    % always report it as running:
    isRunning=1;
else
    % Maybe: Need to scan process-list to report run state of update:
    % N.B.: This seems to be broken on Matlab R2009a and later, but as
    % these releases don't run on ancient OS/X versions anymore, who cares?
    isRunning=length(GetProcessDescriptorFromCommandName('update'));
end
