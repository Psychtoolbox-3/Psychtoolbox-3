function result=EyelinkDoTrackerSetup(el, sendkey)
% USAGE: result=EyelinkDoTrackerSetup(el [, sendkey])
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
% Note that EyelinkDoTrackerSetup() internally uses Beeper() and Snd() to play
% auditory feedback tones if el.targetbeep=1 or el.feedbackbeep=1 and the
% el.callback function is set to the default PsychEyelinkDispatchCallback().
% If you want to use PsychPortAudio in a script that also calls EyelinkDoTrackerSetup,
% then read "help Snd" for instructions on how to provide proper interoperation
% between PsychPortAudio and the feedback sounds created by Eyelink.

%
% 02-06-01  fwc removed use of global el, as suggest by John Palmer.
%               el is now passed as a variable, we also initialize Tracker state bit
%               and Eyelink key values in 'initeyelinkdefaults.m'
% 15-10-02  fwc added sendkey variable that allows to go directly into a particular mode
% 22-06-06  fwc OSX-ed
% 15-06-10  fwc added code for new callback version

if nargin < 1
    error( 'USAGE: result=EyelinkDoTrackerSetup(el [,sendkey])' );
end

% if we have a callback set, we call it.
if ~isempty(el.callback)
    global eyelinkanimationtarget;
    if isempty(eyelinkanimationtarget)
        initmoviestruct();
    end

    global inDoTrackerSetup;
    if isempty(inDoTrackerSetup)
        inDoTrackerSetup = true;
    end

    global inDoDriftCorrection;
    if isempty(inDoDriftCorrection)
        inDoDriftCorrection = false;
    end

    if strcmpi(el.calTargetType, 'video') && ~eyelinkanimationtarget.init
        loadanimationmovie(el);
    end

    result = Eyelink( 'StartSetup', 1 );

    if strcmpi(el.calTargetType, 'video') && ~inDoDriftCorrection && eyelinkanimationtarget.init
        cleanupmovie(el);
    end
    inDoTrackerSetup = false;
    return;
else

    % else we continue with the old version
    if nargin < 2
        sendkey = [];
    end
    result=EyelinkLegacyDoTrackerSetup(el, sendkey);
    return
end


    function initmoviestruct()
        eyelinkanimationtarget.init = false;
        eyelinkanimationtarget.movie = [];
        eyelinkanimationtarget.movieduration = [];
        eyelinkanimationtarget.fps = [];
        eyelinkanimationtarget.imgw = [];
        eyelinkanimationtarget.imgh = [];
        eyelinkanimationtarget.calxy =[];
    end

    function loadanimationmovie(el)
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie', el.window, el.calAnimationTargetFilename, el.calAnimationOpenAsync, el.calAnimationOpenPreloadSecs, el.calAnimationOpenSpecialFlags1);
        eyelinkanimationtarget.init = true;
        eyelinkanimationtarget.movie = movie;
        eyelinkanimationtarget.movieduration = movieduration;
        eyelinkanimationtarget.fps = fps;
        eyelinkanimationtarget.imgw = imgw;
        eyelinkanimationtarget.imgh = imgh;
        Screen('SetMovieTimeIndex', eyelinkanimationtarget.movie, 0, el.calAnimationSetIndexIsFrames);
    end

    function cleanupmovie(el)
        texkill=Screen('GetMovieImage', el.window, eyelinkanimationtarget.movie, el.calAnimationWaitTexClose);
        Screen('PlayMovie', eyelinkanimationtarget.movie, 0, el.calAnimationLoopParam);
        if texkill>0
            Screen('Close', texkill);
        end
        Screen('CloseMovie', eyelinkanimationtarget.movie);
        initmoviestruct();
    end

end

