function [result, dummy]=EyelinkInit(dummy, enableCallbacks)
% USAGE: [result dummy]=EyelinkInit([dummy=0][enableCallbacks=1])
%
% Initialize Eyelink system and connection. Optional arguments:
%
% dummy: Omit, or set to 0 to attempt real initialization,
%        set to 1 to enforce to use initializedummy for dummy mode.
%
%        If regular initialization fails, it will provide option
%        for dummy initilization.
%
% enableCallbacks: Set to 0 for operation without callbacks and
%                  without display of eye camera images on the Subject PC.
%                  Omit or set to 1 for callback and video display operations
%                  during tracker setup, drift correction with the default callback
%                  dispatcher function PsychEyelinkDispatchCallback.m.
%                  Provide namestring of your own dispatcher function if
%                  you want to enable callbacks and video display with a
%                  non-standard, customized dispatcher.
%
% Optional return arguments:
%
% Returns result=1 when succesful, 0 otherwise
% Returns dummy=1 when initialized in dummy mode, 0 otherwise.
%

% history
% 28-10-02	fwc	extracted it from other program
% 09-12-02	fwc	updated to use EyelinkDummyModeDlg instead of Screen('dialog')
%			    a little more testing for erroneous initialization
% 28-06-06  fwc Updated for OSX version
% 21-02-07  fwc Now also returns dummy parameter. Latter has
%               possibly changed.
% 11-04-09  mk  Accept optional 'enableCallbacks' argument. If provided and
%               a string, then use that string as function name for a
%               Eyelink callback function and assign it. If provided and
%               non-zero value, enable callbacks with our default callback
%               function PsychEyelinkDispatchCallback().
%
%               Callbacks are disabled by default to retain backward
%               compatibility.
% 15-06-10 fwc  for consistency, changed PsychEyelinkDispatchCallback to EyelinkDispatchCallback
%               Enabled callbacks by default, because presence of eye image
%               is much desired and expected behaviour
% 15-03-12  mk  Make help text consistent with actual behaviour wrt. to callbacks.

result=0;

% Perform real init by default, if 'dummy' is not provided:
if ~exist('dummy', 'var')
    dummy=[];
end

if isempty(dummy)
    dummy=0;
end

if ~exist('enableCallbacks', 'var')
    % Callbacks enabled by default:
    enableCallbacks = 1;
end

% Callback argument provided and non-empty?
if ~isempty(enableCallbacks)
    % A string with a callback function name? If so, we use it as is:
    if ~ischar(enableCallbacks)
        % Nope. Enable argument 1?
        if enableCallbacks > 0
            % Assign our default callback function:
            enableCallbacks = 'PsychEyelinkDispatchCallback';
        else
            % Callbacks disabled on request:
            enableCallbacks = [];
        end
    end
end

% Dummy init explicitely requested?
if dummy==1
    % Try dummy connection init:
    if Eyelink('InitializeDummy', enableCallbacks) ~=0
        return;
    end
else
    % Try real init of tracker connection:
    if Eyelink('Initialize', enableCallbacks) ~= 0
        % Failed. Ask user if she wants to continue in dummy mode:
        if 1==EyelinkDummyModeDlg
            % Try dummy init:
            if Eyelink('InitializeDummy', enableCallbacks) ~=0
                % Failed as well :-(
                return;
            end

            % Dummy mode active:
            dummy=1;
        else
            % Abort.
            return;
        end
    end
end

result=1;

% Ready:
return;
