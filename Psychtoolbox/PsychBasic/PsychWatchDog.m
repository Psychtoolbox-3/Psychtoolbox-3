function PsychWatchDog(heartbeat)
% PsychWatchDog - Watchdog mechanism and error handler for Psychtoolbox.
%
% PsychWatchDog([heartbeat]);
%
% PsychWatchDog implements a watchdog timer under GNU/Octave. It is not
% supported under the Matlab runtime environment yet.
%
% After PsychWatchDog is enabled, it will become active as soon an Octave
% waits for user input, either via some keyboard input command within a
% script, or when the command prompt becomes active after a script has
% finished running, or been terminated by an error or a users CTRL+C
% interrupt.
%
% It will respond to the user pressing CTRL + . by closing all open
% onscreen windows. If a timeout was set, it will also automatically close
% all windows if the timeout is reached without a sign of life from the
% users script.
%
% Usage:
%
% 1. Call PsychWatchDog(heartbeat) once at the beginning of your script to
% enable the watchdog. A 'heartbeat' value of zero will only enable
% watching for CTRL + . key press. A 'heartbeat' value > 0 will set the
% timeout to 'heartbeat' seconds.
%
% 2. Periodically call PsychWatchDog; without arguments in your script to
% tell the watchdog that your code is still running properly. If Octave
% waits for input and more than 'heartbeat' seconds elapse since the last
% PsychWatchDog call, then the watchdog will trigger and close all windows.
%
% 3. At the end of your script, call PsychWatchDog(-1); to disable the
% watchdog again.
%
% You can call PsychWatchDog() with a vector of at least 2 keycodes as
% delivered by KbName to change the abort key combo from CTRL + . to some
% other arbitrary two key combination. E.g., ...
% PsychWatchDog([KbName('ESCAPE'), KbName('q')]);
% ... would change the abort key combo to ESCape key + q key.
%

% History:
% 1.6.2009   mk Written.

persistent panickeyCodes;
persistent panicdeadline;
persistent heartbeatinterval;

if isempty(panickeyCodes)
    KbName('UnifyKeyNames');    
    panickeyCodes = [ KbName('LeftControl'), KbName('.>') ];
    panicdeadline = inf;
    heartbeatinterval = -1;
    
    % Attach to callback:
    if IsOctave
        add_input_event_hook(mfilename, []);
    end
end

if nargin == 0
    % Heartbeat signal from usercode:
    % Update timeout from heartbeat:
    panicdeadline = GetSecs + heartbeatinterval;
    fprintf('PsychWatchDog: Heartbeat! [Timeout %f] Deadline at %f.\n', heartbeatinterval, panicdeadline);
    return;
end

% Called from usercode to setup/disable/modify settings or heartbeat?
if ~isempty(heartbeat)
    if ~isscalar(heartbeat)
        % heartbeat is a vector of panickeycodes:
        panickeyCodes = heartbeat;
        return;
    end
    
    % Set timeout hearbeat interval value:
    heartbeatinterval = heartbeat;
    panicdeadline     = GetSecs + heartbeatinterval;

    if heartbeat == -1
        % Disable and detach:
        if IsOctave
            remove_input_event_hook(mfilename);
        end
        
        panickeyCodes = [];
    end
    
    return;
end

% Called from runtime system periodically to perform watchdog-cycle:
% Any windows open? Otherwise we got nothing to do...
if (~isempty(Screen('Windows'))) && (heartbeatinterval > -1)
    % Yes. Need to watchdog'em:
    [isdown tnow keycodes] = KbCheck(-1);
    
    if (tnow > panicdeadline) && (heartbeatinterval > 0)
        % Timeout: Close screens.

        % Disable ourselves:
        heartbeatinterval = -1;

        warning('PsychWatchDog: Timeout for heartbeat signal reached. Auto-Closing windows.'); %#ok<WNTAG>
        sca;
        return;
    end
    
    if isdown
        if all(keycodes(panickeyCodes))
            % CTRL + C pressed: Close screens.

            % Disable ourselves:
            heartbeatinterval = -1;

            warning('PsychWatchDog: Abort key combo (e.g., CTRL + .) pressed by user. Auto-Closing windows.'); %#ok<WNTAG>
            sca;
            return;
        end
    end
    
    % Boring watchdog cycle - Everything fine.
    % fprintf('PsychWatchDog: [Timeout %f] All fine at %f.\n', heartbeatinterval, tnow);
end
