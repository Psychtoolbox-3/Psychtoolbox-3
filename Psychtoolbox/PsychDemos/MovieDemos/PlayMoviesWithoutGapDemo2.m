function PlayMoviesWithoutGapDemo2(moviename)
% PlayMoviesWithoutGapDemo2(moviename)
%
% This demo accepts a pattern for a valid moviename, e.g.,
% moviename='*.mpg', then it plays all movies in the current working
% directory whose names match the provided pattern, e.g., the '*.mpg'
% pattern would play all MPEG files in the current directory.
%
% Pressing ESC ends the demo.
%
% Movies are played one after each other. We try to minimize perceptible
% gaps between end of movie i and start of movie i+1 by asynchronous
% loading: While movie i is played back, we ask Psychtoolbox to load the
% next movie i+1 in the background, so startup time for movie i+1 will be
% minimized.
%
% This method uses only one 'movie' object for successive, gapless playback
% of many movie files. The advantage is simplicity in use and maximum
% quality in the sense that movies should play back-to-back without any
% significant delays between them. The downside is loss in flexibility.
% This method only works if all movies have exactly the same image size,
% color depth and format, aspect ratio and playback framerate. They may
% only differ in filename, content and duration. Any other difference will
% cause malfunctions, even hard crashes of Matlab or Octave!
%
% If you need a method that is more flexible and can handle movies of
% different format, use the method in PlayMoviesWithoutGapDemo1. That
% method can handle arbitrary movies. The downside is a higher likelyhood
% of small gaps between movies.
%
% History: 14.03.2012  mk  Wrote it, derived from PlayMoviesWithoutGapDemo1.

AssertOpenGL;

if nargin < 1
    moviename = '*.mov' %#ok<NOPRT>
end

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');
esc=KbName('ESCAPE');

try
    % Open onscreen window:
    screen=max(Screen('Screens'));
    win = Screen('OpenWindow', screen, 128);
    
    Screen('Flip',win);
    abortit = 0;
    
    % Return full list of movie files from directory+pattern:
    moviefiles=dir(moviename);
    for i=1:size(moviefiles,1)
        moviefiles(i).name = [ pwd filesep moviefiles(i).name ];
    end
    
    % Playbackrate defaults to 1:
    rate=1;
    
    % Load first movie. This is a synchronous (blocking) load:
    iteration = 1;
    moviename=moviefiles(mod(iteration, size(moviefiles,1))+1).name;
    [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename);
    
    fprintf('ITER=%i::', iteration);
    fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
    
    % Start playback of movie. This will start
    % the realtime playback clock and playback of audio tracks, if any.
    % Play 'movie', at a playbackrate = rate, with 1.0 == 100% audio volume.
    Screen('PlayMovie', movie, rate, 0, 1.0);
    
    prefetched = 0;
    lastpts = -1;
    
    % Endless loop, runs until ESC key pressed:
    while abortit < 2
        % Show basic info about next movie: Only the name, as we cannot access
        % other info (fps, duration, width and height) for all successor movies.
        % Not a big deal, as all properties have to match the known properties of
        % the first opened movie anyway - except for moviename and duration...
        fprintf('ITER=%i::', iteration);
        fprintf('Movie: %s ...\n', moviename);
        
        i=0;
        
        % Get moviename of next file (after the one that is currently playing):
        iteration=iteration + 1;
        moviename=moviefiles(mod(iteration, size(moviefiles,1))+1).name;
        
        t1 = GetSecs;
        
        % Playback loop: Fetch video frames and display them...
        while 1
            i=i+1;
            if abs(rate) > 0
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex either the texture handle or zero if no new frame is
                % ready yet. pts = Presentation timestamp in seconds.
                [tex pts] = Screen('GetMovieImage', win, movie, 1);
                
                % Valid texture returned?
                if tex < 0
                    % No. This means that the end of this movie is reached.
                    % This can't really happen, unless something went wrong
                    % during playback, because we play all movies round-robin
                    % completely seamless.
                    break
                end
                
                if tex > 0
                    % Yes. Draw the new texture immediately to screen:
                    Screen('DrawTexture', win, tex);
                    
                    % Update display:
                    Screen('Flip', win);
                    
                    % Release texture:
                    Screen('Close', tex);
                end
            end
            
            % Check for abortion by user:
            abortit = 0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown && keyCode(esc))
                % Set the abort-demo flag.
                abortit = 2;
                break;
            end
            
            % We queue the next movie for playback, immediately
            % after start of playback of the current movie, as indicated
            % by the > 0 presentation timestamp:
            if prefetched==0 && pts > 0
                % Queue for background async load operation:
                % We simply set the async flag to 2 and don't query any
                % return values. We pass in the 'movie' handle of the movie
                % which should be succeeded by the new movie 'moviename':
                Screen('OpenMovie', win, moviename, 2, movie);
                prefetched=1;
            end
            
            % Detect when the followup movie has started playback. We detect
            % the change due to a wraparound of the presentation timestamp:
            if prefetched==1 && pts < lastpts
                % New movie has started. Do a new outer-loop iteration to
                % select a new moviefile as successor:
                prefetched = 0;
                lastpts = -1;
                break;
            end
            
            % Keep track of playback time:
            lastpts = pts;
        end
        
        % Print some stats about last played movie:
        telapsed = GetSecs - t1 %#ok<NOPRT,NASGU>
        finalcount=i %#ok<NOPRT,NASGU>
        
        % As playback of the new movie has been started already, we can
        % simply reenter the playback loop:
    end
    
    % End of playback - stop & close the movie:
    Screen('PlayMovie', movie, 0);
    Screen('CloseMovie', movie);
    
    % Close screen:
    Screen('CloseAll');
    
    % Done.
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    sca;
    psychrethrow(psychlasterror);
end
