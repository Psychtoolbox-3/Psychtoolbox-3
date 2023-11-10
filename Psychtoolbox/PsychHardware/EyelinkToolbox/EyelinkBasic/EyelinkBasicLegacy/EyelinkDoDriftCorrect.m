function result=EyelinkDoDriftCorrect(el, x, y, draw, allowsetup)
warning('EyelinkToolbox:LegacyDoDriftCorrect',['The function EyelinkDoDriftCorrect() is deprecated. Please update ', ...
    'your script to use the current method for handling camera setup mode callbacks with PsychEyelinkDispatchCallback.m.']);
warning('off', 'EyelinkToolbox:LegacyDoDriftCorrect');


% USAGE: result=EyelinkLegacyDoDriftCorrect(el [, x, y, draw, allowsetup])
%
% el: eyelink default values
% x,y: position of driftcorrection target
% draw: set to 1 to draw driftcorrection target
% allowsetup: set to 1 to allow to go in to go to trackersetup
%
% Note that EyelinkLegacyDoDriftCorrect() internally uses Beeper() and Snd() to play
% auditory feedback tones if el.targetbeep=1 or el.feedbackbeep=1 and the
% el.callback function is set to the default PsychEyelinkDispatchCallback().
% If you want to use PsychPortAudio in a script that also calls EyelinkDoDriftCorrect,
% then read "help Snd" for instructions on how to provide proper interoperation
% between PsychPortAudio and the feedback sounds created by Eyelink.
%

% /********* PERFORM DRIFT CORRECTION ON TRACKER  *******/
%    /* Performs a drift correction, with target at (x,y). */
%    /* We are explicitly entering a tracker subfunction, */
%    /* so we have to handle link output explicitly. */
%    /* When we finish or abort the drift correction on the tracker, */
%    /* it won't go to another mode until we tell it to. */
%    /* For drift coorection, we can use the */
%    /* drift correction result message to tell when it's done, */
%    /* and what the result was. */

%    /*  Here we display the target ourselves (ignore target updates), */
%    /* wait for local spacebar or for operator trigger or */
%    /* ESC key abort. */
%    /* If operator aborts with ESC, we assume there's a setup */
%    /* problem and go to the setup menu. */

%    /* RETURNS: 0 if OK, 27 if Setup menu was called. */
%

% Eyelink Toolbox version
% 12-05-01    fwc created first version
% 12-05-01    fwc disabled unconditional erasing of screen
% 02-06-01    fwc removed use of global el, as suggested by John Palmer.
% 18-10-02    fwc    made sure missing variables were filled in with defaults
% 15-06-10    fwc added code for new callback version


result=-1; % initialize
if nargin < 1 || ~exist('el', 'var') || isempty(el)
    error( 'USAGE: result=EyelinkLegacyDoDriftCorrect(el [, x, y, draw, allowsetup])' );
end

% fill in missing variables
if ~exist('x', 'var') || ~exist('y', 'var') || isempty(x) || isempty(y)
    [x,y]=WindowCenter(el.window);
end

if ~exist('draw', 'var') || isempty(draw)
    draw=1;
end

if ~exist('allowsetup', 'var') || isempty(allowsetup)
    allowsetup=1;
end

targetrect=[0 0 0 0];

key=1;
while key~= 0
    [key, el]=EyelinkGetKey(el);        % dump old keys
end

if draw==1
    EyelinkLegacyClearCalDisplay(el);        % setup_cal_display()
    targetrect=EyelinkLegacyDrawCalTarget(el, x, y);     % we are told where it should be.
end

if el.targetbeep==1
    EyelinkLegacyCalTargetBeep(el);
end

status=Eyelink( 'DriftCorrStart', x, y);

tickcount=0;
result=el.NO_REPLY;
while result==el.NO_REPLY
    % check for result of drift correction
    result=Eyelink( 'CalResult' );

    [key, el]=EyelinkGetKey(el);        % getkey() HANDLE LOCAL KEY PRESS

    if el.mousetriggersdriftcorr==1 % allow mouse to trigger drift correction (fwc trick)
        [mx,my,button] = GetMouse(el.window);
        if button==1
            if IsInRect(mx,my,targetrect)
                key=el.SPACE_BAR; % fake a key press when mouse is pressed and on target
            end
        end
    end

    switch key
        case el.TERMINATE_KEY,       % breakout key code
            result=el.TERMINATE_KEY;
            return;
        case { 0,  el.JUNK_KEY    }    % No key
        case el.ESC_KEY,   % 27
            if el.allowlocalcontrol==1
                result=el.ESC_KEY;
            end
            if Eyelink('IsConnected') ==-1
                result=el.ESC_KEY;
            end
        case el.SPACE_BAR,                     % 32: we trigger ourselves
            if el.allowlocaltrigger==1
                Eyelink( 'AcceptTrigger');
            end
            if Eyelink('IsConnected') == el.dummyconnected
                result=0;
            end
        otherwise,          % Echo to tracker for remote control
            if el.allowlocalcontrol==1
                Eyelink('SendKeyButton', key, 0, el.KB_PRESS );
            end
    end % switch key
end % while cal_result==NO_REPLY

if draw==1
    EyelinkLegacyEraseCalTarget(el, targetrect); % bit superfluous actually
    EyelinkLegacyClearCalDisplay(el);    % exit_cal_display()
end

if result==el.ESC_KEY || result==-1    % Did we abort drift correction?
    % yes: go to setup menu to fix any problems
    if el.targetbeep==1
        EyelinkLegacyCalDoneBeep(el, 0);
    end
    if allowsetup==1
        EyelinkLegacyDoTrackerSetup(el);
    else
        Eyelink( 'SetOfflineMode');
    end
else
    % Otherwise, we apply the drift correction
    if el.targetbeep==1
        EyelinkLegacyCalDoneBeep(el, 1);
    end
    Eyelink('ApplyDriftCorr' );
    result=0;
end

return;
