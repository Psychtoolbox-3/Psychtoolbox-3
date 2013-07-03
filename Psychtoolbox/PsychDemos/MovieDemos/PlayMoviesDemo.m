function PlayMoviesDemo(moviename, backgroundMaskOut, tolerance, pixelFormat, maxThreads)
% PlayMoviesDemo(moviename [, backgroundMaskOut][, tolerance][, pixelFormat=4][, maxThreads=-1])
%
% This demo accepts a pattern for a valid moviename, e.g.,
% moviename=`*.mpg`, then it plays all movies in the current working
% directory whose names match the provided pattern, e.g., the `*.mpg`
% pattern would play all MPEG files in the current directory.
%
% If you don't specify a moviename, the demo will ask you if it should play
% our standard DualDiscs.mov demo movie, or rather play through a set of
% videos in a playlist which are streamed from the internet. These 'c'ool
% videos may provide you with useful information for your daily work.
%
% This demo uses automatic asynchronous playback for synchronized playback
% of video and sound. Each movie plays until end, then rewinds and plays
% again from the start. Pressing the Cursor-Up/Down key pauses/unpauses the
% movie and increases/decreases playback rate.
% The left- right arrow keys jump in 1 seconds steps. SPACE jumps to the
% next movie in the list. ESC ends the demo.
%
% If the optional RGB color vector backgroundMaskOut is provided, then
% color pixels in the video which are equal or close to backgroundMaskOut will be
% discarded during drawing. E.g., backgroundMaskOut = [255 255 255] would
% discard all white pixels, backgroundMaskOut = [0 0 0] would discard all
% black pixels etc. The optional tolerance parameter allows for some
% lenience, e.g., tolerance = 10 would discard all pixels whose euclidean
% distance in RGB color space is less than 10 units to the backgroundMaskOut
% color. Background color masking requires a graphics card with fragment
% shader support and will fail otherwise.
%
% If the optional `pixelFormat` is specified, it is used to choose
% optimized video playback methods for specific content. Valid values are 1
% or 2 for greyscale video playback, and 7 or 8 for optimized grayscale
% video playback on modern GPU's with GLSL shading support. Values 3, 4, 5
% and 6 play back color video. 4 is the default, 5 or 6 may provide
% significantly improved playback performance on modern GPUs.
%
% If the optional `maxThreads` is specified, it defines the maximum number
% of parallel processing threads that should be used by multi-threaded
% video codecs for playback. A setting of n selects n threads, a setting of
% zero asks to auto-select an optimum number of threads for a given
% computer. By default, a codec specific default number is used, typically
% one thread.
%

% History:
% 10/30/05  mk  Wrote it.
% 07/17/11  mk  Add support for background pixel color removal via shaders.
%               Code cleanup, dead code removal.
% 04/14/12  mk  Code cleanup, refinements for network video streams.
%               Add useful videos to playlist which are less pathetic than
%               the fairy-tales of the iPhone company.
% 06/30/12  mk  Add 'pixelFormat' parameter, disable pathetic Apple PR
%               videos. Enable cool videos on OSX-64Bit as well.
% 08/03/12  mk  Add 'maxThreads' parameter.
% 06/17/13  mk  Add new (c)ool movies, remove Apple PRopaganda videos, cleanup.

theanswer = [];

if (nargin < 1) || isempty(moviename)
    moviename = [];
    theanswer = input('Serious or cool? Type s or c [s/c]? ', 's');
end;

if isempty(moviename)
    moviename = '*.mov';
end

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');

space=KbName('SPACE');
esc=KbName('ESCAPE');
right=KbName('RightArrow');
left=KbName('LeftArrow');
up=KbName('UpArrow');
down=KbName('DownArrow');
shift=KbName('RightShift');

try
    % Child protection
    AssertOpenGL;
    background=[128, 128, 128];
    
    % Open onscreen window:
    screen=max(Screen('Screens'));
    win = Screen('OpenWindow', screen, 0);
    
    shader = [];
    if (nargin > 1) && ~isempty(backgroundMaskOut)
        if nargin < 3
            tolerance = [];
        end
        shader = CreateSinglePassImageProcessingShader(win, 'BackgroundMaskOut', backgroundMaskOut, tolerance);
    end
    
    % Use default pixelFormat if none specified:
    if nargin < 4
        pixelFormat = [];
    end
    
    % Use default maxThreads if none specified:
    if nargin < 5
        maxThreads = [];
    end
    
    % Clear screen to background color:
    Screen('FillRect', win, background);
    
    % Initial display and sync to timestamp:
    Screen('Flip',win);
    iteration = 0;
    abortit = 0;
    
    % Use blocking wait for new frames by default:
    blocking = 1;
    
    % Default preload setting:
    preloadsecs = [];
    
    if isempty(strfind(moviename, 'http'))
        % Return full list of movie files from directory+pattern:
        moviefiles=dir(moviename);
        
        if isempty(moviefiles)
            moviefiles(1).name = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
        else
            for i=1:size(moviefiles,1)
                moviefiles(i).name = [ pwd filesep moviefiles(i).name ];
            end
        end
        
        moviecount = size(moviefiles,1);
    else
        moviefiles(1).name = moviename;
        moviecount = 1;
    end
    
    if strcmpi(theanswer, 'c')
        % Cool stuff, streaming from the web ;-)
        moviefiles = [];
        
        % Make sure a cache directory for buffering exists.
        PsychHomeDir('.cache');
        
        % Increase buffering time for those:
        % No again. This seems to do more harm than good for many streaming sources: preloadsecs = 10;
        
        % The Godfather of Soul giving quick dance lessons:
        moviefiles(end+1).name = 'http://archive.org/download/9lines-JamesBrownDancingLessons151/9lines-JamesBrownDancingLessons151.ogv';
        
        % Randy Pausch's "Last lecture":
        moviefiles(end+1).name = 'http://archive.org/download/LastLecturebyRandyPausch/Last_Lecture_by_Randy_Pausch_Sept_2007_MWV.ogv';
        
        % Richard Stallman talks about the dangers of software patents:
        moviefiles(end+1).name = 'http://archive.org/download/ifso-stallman/ifso-stallman-mpeg1_512kb.mp4';
        
        % Linus Torvalds talks at Aalto University Finnland:
        moviefiles(end+1).name = 'http://archive.org/download/AaltoTalkWithLinusTorvalds/AaltoTalkWithLinusTorvalds.ogv';
        
        % Elon Musk talks about electrical cars, space-flight and solar power:
        moviefiles(end+1).name = 'http://video.ted.com/talk/podcast/2013/None/ElonMusk_2013.mp4';
        moviefiles(end+1).name = 'http://www.oxfordmartin.ox.ac.uk/webcast/201211_musk.mp4';
        
        % FOSDEM 2012 talk about Linux's next generation graphics display server "Wayland":
        moviefiles(end+1).name = 'http://video.fosdem.org/2012/maintracks/k.1.105/Wayland.webm';
        
        % ELC 2012 talk: Gstreamer-1.0 No-longer-compromise-flexibility-for-performance:
        moviefiles(end+1).name = 'http://d17mmld7179ppq.cloudfront.net/gstreamer-10-no-longer-compromise-flexibility-for-performance_52ca47/hd_ready.webm';
        
        % Count all movies in our playlist:
        moviecount = size(moviefiles,2);
        
        if moviecount == 0
            fprintf('Sorry, i do not have any cool movies for your system configuration. Will use boring default movie file.\n');
            moviefiles(1).name = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
        end
        
        % Use polling to wait for new frames when playing movies from the
        % internet. This to make sure we don't time out too early or block
        % for too long if the network connection is slow / high-latency / bad.
        blocking = 0;
        
        % For network playback we use a higher than default caching time:
        % Actually, we don't, as some network streaming videos don't like this, e.g., Youtube: preloadsecs = 10;
    end
    
    % Playbackrate defaults to 1:
    rate=1;
    
    % Choose 16 pixel text size:
    Screen('TextSize', win, 16);
    
    % Endless loop, runs until ESC key pressed:
    while (abortit<2)
        moviename=moviefiles(mod(iteration, moviecount)+1).name;
        iteration=iteration + 1;
        fprintf('ITER=%i::', iteration);
        
        % Show title while movie is loading/prerolling:
        DrawFormattedText(win, moviename, 'center', 'center', [255 255 0], 40);
        Screen('Flip', win);
        
        % Open movie file and retrieve basic info about movie:
        % Some legacy Screen() mex files don't support the pixelFormat or
        % maxThreads flag, therefore we try to do without them if they are not
        % used anyway, avoiding an error abort:
        if ~isempty(pixelFormat) || ~isempty(maxThreads)
            [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename, [], preloadsecs, [], pixelFormat, maxThreads);
        else
            % Legacy compatible call:
            [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename, [], preloadsecs);
        end
        
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, rate, 1, 1.0);
        
        t1 = GetSecs;
        
        % Infinite playback loop: Fetch video frames and display them...
        while 1
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=2;
                break;
            end;
            
            % Check for skip to next movie:
            if (keyIsDown==1 && keyCode(space))
                % Exit while-loop: This will load the next movie...
                break;
            end;
            
            % Only perform video image fetch/drawing if playback is active
            % and the movie actually has a video track (imgw and imgh > 0):
            if ((abs(rate)>0) && (imgw>0) && (imgh>0))
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex is either the positive texture handle or zero if no
                % new frame is ready yet in non-blocking mode (blocking == 0).
                % It is -1 if something went wrong and playback needs to be stopped:
                tex = Screen('GetMovieImage', win, movie, blocking);
                
                % Valid texture returned?
                if tex < 0
                    % No, and there won't be any in the future, due to some
                    % error. Abort playback loop:
                    break;
                end
                
                if tex == 0
                    % No new frame in polling wait (blocking == 0). Just sleep
                    % a bit and then retry.
                    WaitSecs('YieldSecs', 0.005);
                    continue;
                end
                
                % Draw the new texture immediately to screen:
                Screen('DrawTexture', win, tex, [], [], [], [], [], [], shader);
                
                % Update display:
                Screen('Flip', win);
                
                % Release texture:
                Screen('Close', tex);
                
                % Framecounter:
                i=i+1;
            end;
            
            % Further keyboard checks...
            
            if (keyIsDown==1 && keyCode(right))
                % Advance movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) + 1);
            end;
            
            if (keyIsDown==1 && keyCode(left))
                % Rewind movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) - 1);
            end;
            
            if (keyIsDown==1 && keyCode(up))
                % Increase playback rate by 1 unit.
                if (keyCode(shift))
                    rate=rate+0.1;
                else
                    KbReleaseWait;
                    rate=round(rate+1);
                end;
                Screen('PlayMovie', movie, rate, 1, 1.0);
            end;
            
            if (keyIsDown==1 && keyCode(down))
                % Decrease playback rate by 1 unit.
                if (keyCode(shift))
                    rate=rate-0.1;
                else
                    while KbCheck; WaitSecs(0.01); end;
                    rate=round(rate-1);
                end;
                Screen('PlayMovie', movie, rate, 1, 1.0);
            end;
        end;
        
        telapsed = GetSecs - t1;
        fprintf('Elapsed time %f seconds, for %i frames.\n', telapsed, i);
        
        Screen('Flip', win);
        KbReleaseWait;
        
        % Done. Stop playback:
        Screen('PlayMovie', movie, 0);
        
        % Close movie object:
        Screen('CloseMovie', movie);
    end;
    
    % Close screens.
    Screen('CloseAll');
    
    % Done.
    return;
catch %#ok<*CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
end
