function result=EyelinkLegacyDoTrackerSetup(el, sendkey)
warning('EyelinkToolbox:LegacyDoTrackerSetup',['Use of the function EyelinkDoTrackerSetup() without providing a callback handler ', ...
    '(such as the included PsychEyelinkDispatchCallback) is deprecated. Please update your script to use the currently supported conventions.']);
warning('off', 'EyelinkToolbox:LegacyDoTrackerSetup');

% USAGE: result=EyelinkLegacyDoTrackerSetup(el [, sendkey])
%
% el: Eyelink default values
%
% sendkey:  set to go directly into a particular mode
%           sendkey is optional and ignored if el.callback is defined for
%           callback based tracker setup.
%
%           'v', start validation
%           'c', start calibration
%           'd', start driftcorrection
%           13, or el.ENTER_KEY, show 'eye' setup image
%
% Note that EyelinkLegacyDoTrackerSetup() internally uses Beeper() and Snd() to play
% auditory feedback tones if el.targetbeep=1 or el.feedbackbeep=1 and the
% el.callback function is set to the default PsychEyelinkDispatchCallback().
% If you want to use PsychPortAudio in a script that also calls EyelinkLegacyDoTrackerSetup,
% then read "help Snd" for instructions on how to provide proper interoperation
% between PsychPortAudio and the feedback sounds created by Eyelink.

%
% 02-06-01  fwc removed use of global el, as suggest by John Palmer.
%               el is now passed as a variable, we also initialize Tracker state bit
%               and Eyelink key values in 'initeyelinkdefaults.m'
% 15-10-02  fwc added sendkey variable that allows to go directly into a particular mode
% 22-06-06  fwc OSX-ed
% 15-06-10  fwc added code for new callback version

result=-1;
if nargin < 1
    error( 'USAGE: result=EyelinkLegacyDoTrackerSetup(el [,sendkey])' );
end

% if we have the new callback code, we call it.
if ~isempty(el.callback)
    error('el.callback is not empty. Legacy functions not supported when callback is set.');
end

Eyelink( 'StartSetup' );        % start setup mode
Eyelink( 'WaitForModeReady', el.waitformodereadytime );  % time for mode change

EyelinkLegacyClearCalDisplay(el);    % setup_cal_display()
key=1;
while key~= 0
    key=EyelinkGetKey(el);        % dump old keys
end

% go directly into a particular mode
if nargin==2 && ~isempty(sendkey)
    if el.allowlocalcontrol==1
        switch lower(sendkey)
            case{ 'c', 'v', 'd', el.ENTER_KEY}
                forcedkey=double(sendkey(1,1));
                Eyelink('SendKeyButton', forcedkey, 0, el.KB_PRESS );
        end
    end
end

tstart=GetSecs;
stop=0;
while stop==0 && bitand(Eyelink( 'CurrentMode'), el.IN_SETUP_MODE)

    i=Eyelink( 'CurrentMode');

    if ~Eyelink( 'IsConnected' )
        stop=1;
        break;
    end

    if bitand(i, el.IN_TARGET_MODE)            % calibrate, validate, etc: show targets
        EyelinkLegacyTargetModeDisplay(el);
    elseif bitand(i, el.IN_IMAGE_MODE)        % display image until we're back
        if Eyelink ('ImageModeDisplay')==el.TERMINATE_KEY
            result=el.TERMINATE_KEY;
            return;    % breakout key pressed
        else
            EyelinkLegacyClearCalDisplay(el);    % setup_cal_display()
        end
    end

    [key, el]=EyelinkGetKey(el);        % getkey() HANDLE LOCAL KEY PRESS
    if 1 && key~=0 && key~=el.JUNK_KEY    % print pressed key codes and chars
        fprintf('%d\t%s\n', key, char(key) );
    end

    switch key
        case el.TERMINATE_KEY                % breakout key code
            result=el.TERMINATE_KEY;
            return;
        case { 0, el.JUNK_KEY }          % No or uninterpretable key
        case el.ESC_KEY
            if Eyelink('IsConnected') == el.dummyconnected
                stop=1; % instead of 'goto exit'
            end
            if el.allowlocalcontrol==1
                Eyelink('SendKeyButton', key, 0, el.KB_PRESS );
            end
        otherwise         % Echo to tracker for remote control
            if el.allowlocalcontrol==1
                Eyelink('SendKeyButton', double(key), 0, el.KB_PRESS );
            end
    end
end % while IN_SETUP_MODE

% exit:
EyelinkLegacyClearCalDisplay(el);    % exit_cal_display()
result=0;
return;
