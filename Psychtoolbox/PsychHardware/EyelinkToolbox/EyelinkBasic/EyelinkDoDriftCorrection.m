function success=EyelinkDoDriftCorrection(el, x, y, draw, allowsetup)
% success=EyelinkDoDriftCorrection(el [, x][, y][, draw][, allowsetup])
%
% DO PRE-TRIAL DRIFT CORRECTION
% We repeat if ESC key pressed to do setup.
% Setup might also have erased any pre-drawn graphics.
%
% Note that EyelinkDoDriftCorrection() internally uses Beeper() and Snd() to play
% auditory feedback tones if el.targetbeep=1 or el.feedbackbeep=1 and the
% el.callback function is set to the default PsychEyelinkDispatchCallback().
% If you want to use PsychPortAudio in a script that also calls EyelinkDoDriftCorrection,
% then read "help Snd" for instructions on how to provide proper interoperation
% between PsychPortAudio and the feedback sounds created by Eyelink.
%

success=1;

% if no x and y are supplied, set x,y to center coordinates
if ~exist('x', 'var') || isempty(x) || ~exist('y', 'var') || isempty(y)
    [x,y] = WindowCenter(el.window); % convenience routine part of eyelink toolbox
end

if ~exist('draw', 'var') || isempty(draw)
    draw=1;
end

if ~exist('allowsetup', 'var') || isempty(allowsetup)
    allowsetup=1;
end

while 1
    if Eyelink('IsConnected')==el.notconnected   % Check link often so we don't lock up if tracker lost
        success=0;
        return;
    end
    
    if ~isempty(el.callback) % if we have a callback set, we call it.
        global eyelinkanimationtarget;
        if isempty(eyelinkanimationtarget)
            initmoviestruct();
        end

        global inDoTrackerSetup;
        if isempty(inDoTrackerSetup)
            inDoTrackerSetup = false;
        end

        global inDoDriftCorrection;
        if isempty(inDoDriftCorrection)
            inDoDriftCorrection = true;
        end

        if strcmpi(el.calTargetType, 'video') && ~eyelinkanimationtarget.init
            loadanimationmovie(el);
        end
        result = Eyelink('DriftCorrStart', x, y, 1, draw, allowsetup);

    else
        % else we continue with the old version
        result = EyelinkLegacyDoDriftCorrect(el, x, y, draw, allowsetup);
    end
    
    if result==el.TERMINATE_KEY
        success=0;
        return;
    end
    
    % repeat if ESC was pressed to access Setup menu
    if(result~=el.ESC_KEY)
        break;
    end
    
end % while

% fprintf('~isempty(el.callback): %d, el.calTargetType: %s, ~inDoTrackerSetup: %d, eyelinkanimationtarget.init: %d\n', ~isempty(el.callback), el.calTargetType, ~inDoTrackerSetup, eyelinkanimationtarget.init)
if ~isempty(el.callback) && strcmpi(el.calTargetType, 'video') && ~inDoTrackerSetup && eyelinkanimationtarget.init
    cleanupmovie(el);
    inDoDriftCorrection = false;
end

return

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
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie',  el.window, el.calAnimationTargetFilename, el.calAnimationOpenAsync, el.calAnimationOpenPreloadSecs, el.calAnimationOpenSpecialFlags1);
        eyelinkanimationtarget.init = true;
        eyelinkanimationtarget.movie = movie;
        eyelinkanimationtarget.movieduration = movieduration;
        eyelinkanimationtarget.fps = fps;
        eyelinkanimationtarget.imgw = imgw;
        eyelinkanimationtarget.imgh = imgh;
        Screen('SetMovieTimeIndex', eyelinkanimationtarget.movie, 0, el.calAnimationSetIndexIsFrames);
    end

    function cleanupmovie(el)
        texkill = Screen('GetMovieImage', el.window, eyelinkanimationtarget.movie, el.calAnimationWaitTexClose);
        Screen('PlayMovie', eyelinkanimationtarget.movie, 0, el.calAnimationLoopParam);
        if texkill > 0
            Screen('Close', texkill);
        end
        Screen('CloseMovie', eyelinkanimationtarget.movie);
        initmoviestruct();
    end
end
