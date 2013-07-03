function PlayMoviesWithoutGapDemo1(moviename)
%
% PlayMoviesWithoutGapDemo1(moviename)
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

% History:
% 7/5/06  mk  Wrote it.

AssertOpenGL;

if nargin < 1
    moviename = '*.mov' %#ok<NOPRT>
end;

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');

esc=KbName('ESCAPE');

try
    % Child protection
    AssertOpenGL;
    background=[128, 128, 128];
    
    % Open onscreen window:
    screen=max(Screen('Screens'));
    win = Screen('OpenWindow', screen, background);
    
    % Initial display and sync to timestamp:
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
    
    % Start playback of movie. This will start
    % the realtime playback clock and playback of audio tracks, if any.
    % Play 'movie', at a playbackrate = 1, with 1.0 == 100% audio volume.
    Screen('PlayMovie', movie, rate, 0, 1.0);
    
    prefetched=0;
    newmovie = -1;
    
    % Endless loop, runs until ESC key pressed:
    while (abortit<2)
        % Show basic info about movie:
        fprintf('ITER=%i::', iteration);
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        
        % Get moviename of next file:
        iteration=iteration + 1;
        moviename=moviefiles(mod(iteration, size(moviefiles,1))+1).name;
        
        t1 = GetSecs;
        
        % Playback loop: Fetch video frames and display them...
        while(1)
            i=i+1;
            if (abs(rate)>0)
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex either the texture handle or zero if no new frame is
                % ready yet. pts = Presentation timestamp in seconds.
                [tex pts] = Screen('GetMovieImage', win, movie, 1);
                
                % Valid texture returned?
                if tex<0
                    % No. This means that the end of this movie is reached.
                    % We exit the loop and prepare next movie:
                    break;
                end;
                
                if (tex>0)
                    % Yes. Draw the new texture immediately to screen:
                    Screen('DrawTexture', win, tex);
                    
                    % Update display:
                    Screen('Flip', win);
                    
                    % Release texture:
                    Screen('Close', tex);
                end;
            end;
            
            % Check for abortion by user:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=2;
                break;
            end;
            
            % We start background loading of the next movie 0.5 seconds
            % after start of playback of the current movie:
            if prefetched==0 && pts > 0.5
                % Initiate background async load operation:
                % We simply set the async flag to 1 and don't query any
                % return values:
                Screen('OpenMovie', win, moviename, 1);
                prefetched=1;
            end;
            
            % If asynchronous load of next movie has been started already
            % and we are less than 0.5 seconds from the end of the current
            % movie, then we try to finish the async load operation and
            % start playback of the new movie in order to give processing a
            % headstart:
            if prefetched==1 && movieduration - pts < 0.5
                % Less than 0.5 seconds until end of current movie. Try to
                % start playback for next movie:
                
                % Open the movie - this should take zero-time as the movie
                % should have been already prefetched during playback of
                % the current movie. In case loading did not finish, this
                % command will pause Matlabs execution until the movie is
                % really ready.
                [newmovie newmovieduration newfps newimgw newimgh] = Screen('OpenMovie', win, moviename);
                
                % Start it:
                Screen('PlayMovie', newmovie, rate, 0, 1.0);
                
                prefetched=2;
            end;
        end;
        
        telapsed = GetSecs - t1 %#ok<NOPRT,NASGU>
        finalcount=i %#ok<NOPRT,NASGU>
        
        % Done with old movie. Stop its playback:
        Screen('PlayMovie', movie, 0);
        
        % Close movie object:
        Screen('CloseMovie', movie);
        
        % Ok, now our 'newmovie' becomes the current 'movie':
        movie = newmovie;
        movieduration = newmovieduration;
        fps = newfps;
        imgw = newimgw;
        imgh = newimgh;
        prefetched = 0;
        % As playback of the new movie has been started already, we can
        % simply reenter the playback loop:
    end;
    
    % End of movie playback. Shut down and exit:
    if prefetched == 1
        % A prefetch operation for a movie is still in progress. We need to
        % finalize this cleanly by waiting for the movie to open and then
        % closing it.
        newmovie = Screen('OpenMovie', win, moviename);
        prefetched = 2;
    end
    
    if prefetched == 2
        % Playback of a new prefetched movie has been started. We need to stop and
        % close it:
        Screen('PlayMovie', newmovie, 0);
        Screen('CloseMovie', newmovie);
    end
    
    % Close screens.
    Screen('Close', win);
    
    % Done.
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
end;
