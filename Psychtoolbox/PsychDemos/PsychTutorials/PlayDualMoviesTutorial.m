function PlayDualMoviesTutorial(moviename)
%
% PlayDualMoviesTutorial(moviename)
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
% 06/17/13  mk  Cleaned up.

if nargin < 1
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

try
    % Child protection
    AssertOpenGL;
    
    % Open onscreen window:
    screen=max(Screen('Screens'));
    [win, scr_rect] = Screen('OpenWindow', screen, 0, [], [], [], 0, [], kPsychNeedFastBackingStore);
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % We create a Luminance+Alpha matrix for use as transparency mask:
    % Layer 1 (Luminance) is filled with 'backgroundcolor'.
    transLayer=2;
    ms=200;
    [x,y]=meshgrid(-ms:ms, -ms:ms);
    maskblob=ones(2*ms+1, 2*ms+1, transLayer) * 255;
    % Layer 2 (Transparency aka Alpha) is filled with gaussian transparency
    % mask.
    xsd=ms/2.2;
    ysd=ms/2.2;
    maskblob(:,:,transLayer)=round(255 - exp(-((x/xsd).^2)-((y/ysd).^2))*255);
    
    % Build a single transparency mask texture:
    masktex=Screen('MakeTexture', win, maskblob);
    
    % Initial display and sync to timestamp:
    Screen('Flip',win);
    iteration=2;
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
        [xc yc] = RectCenter(rect2);
        SetMouse(xc, yc, win);
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
        
        angle=0;
        alpha=1;
        
        % Infinite playback loop: Fetch video frames and display them...
        while 1
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
                    Screen('SelectStereoDrawBuffer', win, 0);
                    Screen('DrawTexture', win, tex, [], rect1);
                    % Release texture:
                    Screen('Close', tex);
                end;
                
                % Valid 2nd texture returned?
                if tex2>0
                    % Draw the new texture immediately to screen:
                    Screen('SelectStereoDrawBuffer', win, 1);
                    Screen('DrawTexture', win, tex2, [], rect2, angle, [], alpha);
                    % Release texture:
                    Screen('Close', tex2);
                end;
                
                % Update display if there is anything to update:
                if (tex>0 || tex2>0)
                    % We use clearmode=1, aka don't clear on flip. This is
                    % needed to avoid flicker...
                    Screen('Flip', win, 0, 0);
                end;
            end;
            
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown==1 && keyCode(esc))
                break;
            end;
            
            if (keyIsDown==1 && keyCode(space))
                % Exit while-loop: This will load the next movie...
                KbReleaseWait;
                break;
            end;
            
            if (keyIsDown==1 && keyCode(right))
                angle = angle + 1;
            end;
            
            if (keyIsDown==1 && keyCode(left))
                angle = angle - 1;
            end;
            
            if (keyIsDown==1 && keyCode(up))
                alpha = alpha + 0.01;
                if alpha > 1
                    alpha = 1;
                end
            end;
            
            if (keyIsDown==1 && keyCode(down))
                alpha = alpha - 0.01;
                if alpha < 0
                    alpha = 0;
                end
            end;
            
            [x y] = GetMouse(win);
            rect2 = CenterRectOnPoint(rect2, x, y);
        end;
        
        rect2 = rect1;
        
        % Infinite playback loop: Fetch video frames and display them...
        while 1
            i=i+1;
            if (abs(rate)>0)
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex either the texture handle or zero if no new frame is
                % ready yet.
                tex = Screen('GetMovieImage', win, movie, 1);
                tex2 = Screen('GetMovieImage', win, movie2, 1);
                
                Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                Screen('DrawTexture', win, masktex, [], rect2);
                
                % Valid texture returned?
                if tex>0
                    % Draw the new texture immediately to screen:
                    Screen('BlendFunction', win, GL_DST_ALPHA, GL_ZERO);
                    Screen('DrawTexture', win, tex, [], rect1);
                    % Release texture:
                    Screen('Close', tex);
                end;
                
                % Valid 2nd texture returned?
                if tex2>0
                    % Overdraw the first tex with 2nd tex:
                    Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_ONE);
                    Screen('DrawTexture', win, tex2, rect2, rect2);
                    % Release texture:
                    Screen('Close', tex2);
                end;
                
                % Update display if there is anything to update:
                if (tex>0 || tex2>0)
                    % We use clearmode=1, aka don't clear on flip. This is
                    % needed to avoid flicker...
                    Screen('Flip', win, 0, 0);
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
                angle = angle + 1;
            end;
            
            if (keyIsDown==1 && keyCode(left))
                angle = angle - 1;
            end;
            
            if (keyIsDown==1 && keyCode(up))
                alpha = alpha + 0.01;
                if alpha > 1
                    alpha = 1;
                end
            end;
            
            if (keyIsDown==1 && keyCode(down))
                alpha = alpha - 0.01;
                if alpha < 0
                    alpha = 0;
                end
            end;
            
            [x y] = GetMouse(win);
            rect2 = ClipRect(CenterRectOnPoint([0 0 401 401], x, y), scr_rect);
        end;
        
        Screen('Flip', win);
        Screen('Flip', win);
        
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
end;
