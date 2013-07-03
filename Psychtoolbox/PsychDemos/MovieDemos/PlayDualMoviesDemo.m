function PlayDualMoviesDemo(moviename)
%
% PlayDualMoviesDemo(moviename)
%
% This demo accepts a pattern for a valid moviename, e.g.,
% moviename='*.mpg', then it plays all movies in the current working
% directory whose names match the provided pattern, e.g., the '*.mpg'
% pattern would play all MPEG files in the current directory.
%
% This demo uses automatic asynchronous playback for synchronized playback
% of video and sound. Each movie plays until end, then rewinds and plays
% again from the start. Pressing the Cursor-Up/Down key pauses/unpauses the
% movie and increases/decreases playback rate.
% The left- right arrow keys jump in 1 seconds steps. SPACE jumps to the
% next movie in the list. ESC ends the demo.
%

% History:
% 10/30/05  mk  Wrote it.
% 06/17/13  mk  Cleanup.

if nargin < 1
    moviename = '*.mov';
end;

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
    
    % Grey background:
    background=[128, 128, 128];
    
    % Open onscreen window:
    screen=max(Screen('Screens'));
    win = Screen('OpenWindow', screen, background);
    
    % Initial display and sync retrace:
    Screen('Flip',win);
    iteration=0;
    abortit=0;
    
    % Return full list of movie files from directory+pattern:
    moviefiles=dir(moviename);
    for i=1:size(moviefiles,1)
        moviefiles(i).name = [ pwd filesep moviefiles(i).name ];
    end
    
    % Endless loop, runs until ESC key pressed:
    while (abortit<2)
        iteration=iteration + 1;
        moviename=moviefiles(mod(iteration, size(moviefiles,1))+1).name;
        moviename2=moviefiles(mod(iteration+1, size(moviefiles,1))+1).name;
        
        % Open movie file and retrieve basic info about movie:
        [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename);
        fprintf('Movie1: %s  : %f seconds duration, %f fps...\n', moviename, movieduration, fps);
        rect1=SetRect(1,1,imgw,imgh);
        % Open 2nd movie file and retrieve basic info about movie:
        [movie2 movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename2);
        fprintf('Movie2: %s  : %f seconds duration, %f fps...\n', moviename2, movieduration, fps);
        rect2=SetRect(1,1,imgw,imgh);
        rect2=AdjoinRect(rect2, rect1, RectRight);
        
        i=0;
        
        % Seek to start of movies (timeindex 0):
        Screen('SetMovieTimeIndex', movie, 0);
        Screen('SetMovieTimeIndex', movie2, 0);
        
        rate=1;
        
        % Start playback of movies. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, rate, 1, 1.0);
        Screen('PlayMovie', movie2, rate, 1, 1.0);
        
        t1 = GetSecs;
        
        % Infinite playback loop: Fetch video frames and display them...
        while(1)
            i=i+1;
            if (abs(rate)>0)
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex either the texture handle or zero if no new frame is
                % ready yet.
                tex = Screen('GetMovieImage', win, movie, 1);
                tex2 = Screen('GetMovieImage', win, movie2, 1);
                
                % Valid texture returned?
                if tex>0
                    % Draw the new texture immediately to screen:
                    Screen('DrawTexture', win, tex, [], rect1);
                    % Release texture:
                    Screen('Close', tex);
                end;
                
                % Valid 2nd texture returned?
                if tex2>0
                    % Draw the new texture immediately to screen:
                    Screen('DrawTexture', win, tex2, [], rect2);
                    % Release texture:
                    Screen('Close', tex2);
                end;
                
                % Update display if there is anything to update:
                if (tex>0 || tex2>0)
                    % We use clearmode=1, aka don't clear on flip. This is
                    % needed to avoid flicker...
                    Screen('Flip', win, 0, 1);
                end;
            end;
            
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=2;
                break;
            end;
            
            if (keyIsDown==1 && keyCode(space))
                % Exit while-loop: This will load the next movie...
                break;
            end;
            
            if (keyIsDown==1 && keyCode(right))
                % Advance movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) + 1);
                Screen('SetMovieTimeIndex', movie2, Screen('GetMovieTimeIndex', movie2) + 1);
            end;
            
            if (keyIsDown==1 && keyCode(left))
                % Rewind movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) - 1);
                Screen('SetMovieTimeIndex', movie2, Screen('GetMovieTimeIndex', movie2) - 1);
            end;
            
            if (keyIsDown==1 && keyCode(up))
                % Increase playback rate by 1 unit.
                KbReleaseWait;
                if (keyCode(shift))
                    rate=rate+0.1;
                else
                    rate=round(rate+1);
                end;
                Screen('PlayMovie', movie, rate, 1, 1.0);
                Screen('PlayMovie', movie2, rate, 1, 1.0);
            end;
            
            if (keyIsDown==1 && keyCode(down))
                % Decrease playback rate by 1 unit.
                KbReleaseWait;
                if (keyCode(shift))
                    rate=rate-0.1;
                else
                    rate=round(rate-1);
                end;
                Screen('PlayMovie', movie, rate, 1, 1.0);
                Screen('PlayMovie', movie2, rate, 1, 1.0);
            end;
        end;
        
        fprintf('Elapsed time %f secs, count %i.\n', GetSecs - t1, i);
        
        Screen('Flip', win);
        Screen('Flip', win);
        
        % Wait for all keys released:
        KbReleaseWait;
        
        % Done. Stop playback:
        Screen('PlayMovie', movie, 0);
        Screen('PlayMovie', movie2, 0);
        
        % Close movie objects:
        Screen('CloseMovie', movie);
        Screen('CloseMovie', movie2);
    end;
    
    % Close screens.
    Screen('CloseAll');
    
    % Done.
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
end
