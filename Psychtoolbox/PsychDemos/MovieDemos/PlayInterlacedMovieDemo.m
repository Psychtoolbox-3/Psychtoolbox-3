function PlayInterlacedMovieDemo(moviename)
%
% PlayInterlacedMovieDemo(moviename)
%
% A simple demo to show playback of interlaced movies. Input
% images from the movie are deinterlaced on-the-fly by use of a simple GLSL
% deinterlace shader.
%
% The demo reads a movie file from the file 'moviename' and plays it once,
% then exits.
%
% Should work with NVidia GeforceFX 5200 or ATI Radeon 9600 and later, or
% Intel GMA 950.

% History:
% 10/30/05  mk  Wrote it.
% 06/17/13  mk  Cleaned up.

if nargin < 1
    error('Moviename missing! You need to provide the filename of a movie file.');
end;

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');
esc=KbName('ESCAPE');
space=KbName('space');

try
    % Child protection
    AssertOpenGL;
    
    % Open onscreen window with black background:
    screen=max(Screen('Screens'));
    win = Screen('OpenWindow', screen, 0);
    
    % Load deinterlacing-shader: This will abort our script if the graphics
    % hardware doesn't support GLSL or fragment shaders.
    deinterlacer = LoadGLSLProgramFromFiles('EXPDeinterlaceShaderLineDouble',1);
    
    % Bind and initialize deinterlacer:
    glUseProgram(deinterlacer);
    
    % Input image will be bound to texture unit zero:
    glUniform1i(glGetUniformLocation(deinterlacer, 'Image1'), 0);
    
    % Get handle for the field selection parameter:
    useoddfield=glGetUniformLocation(deinterlacer, 'UseOddField');
    
    
    % Initial display and sync to retrace:
    Screen('Flip',win);
    iteration=0;
    abortit=0;
    switchfielddisplayorder=0;
    
    % Endless loop, runs until ESC key pressed.
    while ~abortit
        iteration=iteration + 1;
        
        % Open movie file and retrieve basic info about movie:
        [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename);
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, without looping = 0 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, 1, 0, 1.0);
        
        % Infinite playback loop: Fetch video frames and display them...
        while 1
            i=i+1;
            
            % Return next frame in movie, in sync with current playback
            % time and sound.
            % tex is either the texture handle or zero if no new frame is
            % ready yet, or -1 if the movie finished playback.
            % pts = Presentation timestamp in seconds.
            tex = Screen('GetMovieImage', win, movie, 1);
            
            % Valid texture returned?
            if tex<=0
                % Movie has reached its end. Exit our playback loop.
                break;
            end
            
            % Select even half-field (the even rows) in the interlaced
            % texture: One can switch the order of drawing by toggling
            % switchfielddisplayorder between 0 and 1 via press of the
            % space key.
            glUseProgram(deinterlacer);
            glUniform1f(useoddfield, switchfielddisplayorder);
            
            % Draw the new texture immediately to screen, appying the
            % deinterlacer. You need to disable bilinear filtering via the
            % 0 - flag, because the deinterlacer doesn't work yet with
            % anything else than nearest neighbour filtering:
            Screen('DrawTexture', win, tex, [], [], [], 0, [], [], deinterlacer);
            
            % Update display at next monitor retrace:
            Screen('Flip', win);
            
            % Select odd half-field (the odd rows) in the interlaced
            % texture:
            glUseProgram(deinterlacer);
            glUniform1f(useoddfield, 1-switchfielddisplayorder);
            
            % Draw the new texture immediately to screen:
            Screen('DrawTexture', win, tex, [], [], [], 0, [], [], deinterlacer);
            
            % Update display:
            Screen('Flip', win);
            
            % Release texture:
            Screen('Close', tex);
            
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck; %#ok<ASGLU>
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=1;
                break;
            end
            
            if (keyIsDown==1 && keyCode(space))
                % Toggle display order of even and odd fields:
                switchfielddisplayorder = 1 - switchfielddisplayorder %#ok<NOPRT>
                KbReleaseWait;
            end
        end
        
        Screen('Flip', win);
        
        % Done. Stop playback:
        Screen('PlayMovie', movie, 0);
        
        % Close movie object:
        Screen('CloseMovie', movie);
        
        % ...and repeat...
    end
    
    % Close screens.
    Screen('CloseAll');
    
    % Done.
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
end
