function DetectionRTInVideoDemo(moviename, timeOfEvent, trials)
%
% DetectionRTInVideoDemo(moviename, timeOfEvent, trials)
%
% A demo implementation of how to collect reaction time for detection of a
% specific time-locked event in a movie file.
%
% Parameters:
% moviename - Filename of moviefile to use. If none is provided, then the
% simple DualDiscs collision movie is used which is part of PTB.
%
% timeOfEvent - Time (in seconds) when the timelocked event happens which
% should be detected by the subject. Default is time of contact of the
% DualDiscs collision.
%
% trials - Number of trials to run. Defaults to 10 or ESC key press.
%
% How the demo works: Read the source code - its well documented ;-)
% This demo demonstrates the optimal way if you have movies with sound
% (that needs to be played in sync with the video) or long movies that
% don't fit into memory or that you don't want to load into memory.
%
% It uses automatic synchronized audio-video playback. The central
% while-loop checks for arrival of new frames for display and displays them
% when their time has come. It also registers subjects keypress responses
% and does all the bookkeeping and sanity checking.
%
% An alternative approach would be to preload the whole movie into textures:
% The whole movie gets read into PTB textures before start of trial. Then
% you show the textures in quick succession like in MovieDemo. The
% advantage of that approach is exact control over display timing and
% display order: You can show frames in any order you like at any rate you
% like (and that your hardware likes). Disadvantage would be: Longer trial
% setup time for loading the whole movie, higher memory consumption (keep
% all n frames in memory instead of only the current one) and the inability
% to play sound in sync with video.
%

% History:
% 12/19/05  mk  Wrote it.
% 02/03/06  mk  Adapted for use on Windows.
% 03/11/12  mk  Cleanup.
% 06/17/13  mk  Cleanup.

if nargin < 1
    % Default movie is our own disc collision movie:
    moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
end;

if nargin < 2
    timeOfEvent = 1.16;
end;

if nargin < 3
    trials = 10;
end;

fprintf('Loading probe movie %s ...\n', moviename);
fprintf('Time locked probe event will happen at time %f secs after start of movie.\n', timeOfEvent);
fprintf('Will run %i trials or until ESCape key is pressed.\n', trials);

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');

% Query keycodes for ESCAPE key and Space key:
esc=KbName('ESCAPE');
space=KbName('space');

try
    % Child protection: Make sure we run on the OSX / OpenGL Psychtoolbox.
    % Abort if we don't:
    AssertOpenGL;
    
    % Background color will be a grey one:
    background=[128, 128, 128];
    
    % Open onscreen window. We use the display with the highest number on
    % multi-display setups:
    screen=max(Screen('Screens'));
    
    % This will open a screen with background color 'background':
    win = Screen('OpenWindow', screen, background);
    
    % Hide the mouse cursor:
    HideCursor;
    
    % Show instructions...
    tsize=30;
    Screen('TextSize', win, tsize);
    [x, y]=Screen('DrawText', win, 'Collision detection fake experiment.',40, 100); %#ok<*ASGLU>
    [x, y]=Screen('DrawText', win, 'Press ESC-ape key to abort anytime.', 40, y + 10 + tsize);
    [x, y]=Screen('DrawText', win, 'Press SPACE key when you see the discs colliding', 40, y + 10 + tsize);
    Screen('DrawText', win, 'Press any key to start the experiment...', 40, y + 10 + tsize);
    
    % Flip to show the grey screen:
    Screen('Flip',win);
    
    % Wait for keypress + release...
    KbStrokeWait;
    
    % Show cleared screen...
    Screen('Flip',win);
    
    % Wait a second...
    WaitSecs(1);
    
    % Main trial loop: Do 'trials' trials...
    for i=1:trials
        % Open the moviefile and query some infos like duration, framerate,
        % width and height of video frames. We could also query the total count of frames in
        % the movie, but computing 'framecount' takes long, so avoid to query
        % this property if you don't need it!
        [movie movieduration fps] = Screen('OpenMovie', win, moviename);
        
        % We estimate framecount instead of querying it - faster:
        framecount = movieduration * fps;
        
        % Start playback of the movie:
        % Play 'movie', at a playbackrate = 1 (normal speed forward),
        % play it once, aka with loopflag = 0,
        % play audio track at volume 1.0  = 100% audio volume.
        Screen('PlayMovie', movie, 1, 0, 1.0);
        
        % Video playback and key response RT collection loop:
        % This loop repeats until either the subject responded with a
        % keypress to indicate s(he) detected the event in the vido, or
        % until the end of the movie is reached.
        movietexture=0;     % Texture handle for the current movie frame.
        reactiontime=-1;    % Variable to store reaction time.
        lastpts=0;          % Presentation timestamp of last frame.
        onsettime=-1;       % Realtime at which the event was shown to the subject.
        rejecttrial=0;      % Flag which is set to 1 to reject an invalid trial.
        
        while(movietexture>=0 && reactiontime==-1)
            % Check if a new movie video frame is ready for visual
            % presentation: This call polls for arrival of a new frame. If
            % a new frame is ready, it converts the video frame into a
            % Psychtoolbox texture image and returns a handle in
            % 'movietexture'. 'pts' contains a so called presentation
            % timestamp. That is the time (in seconds since start of movie)
            % at which this video frame should be shown on the screen.
            % Arrival of textures is automatically synchronized to the
            % audio track and to real-world time. If the video display loop
            % can't keep up with the flow of time and the soundtrack,
            % the engine will automatically skip/drop frames to keep video
            % in sync with audio as good as possible. If the pts of a new
            % texture is greater than the 'timeOfEvent' then you'll know
            % that this texture will show the visual target event as soon
            % as you draw and 'Flip' it.
            % In case that no new video texture is available yet for
            % presentation, this function will return a zero texture handle
            % to indicate this. If no new texture will become available
            % anymore, because the end of the movie is reached, it will
            % return a handle of -1 to indicate end of playback.
            
            % The 0 - flag means: Don't wait for arrival of new frame, just
            % return a zero or -1 'movietexture' if none is ready.
            [movietexture pts] = Screen('GetMovieImage', win, movie, 0);
            
            % Is it a valid texture?
            if (movietexture>0)
                % Yes. Draw the texture into backbuffer:
                Screen('DrawTexture', win, movietexture);
                
                % Flip the display to show the image at next retrace:
                % vbl will contain the exact system time of image onset on
                % screen: This should be accurate in the sub-millisecond
                % range.
                vbl=Screen('Flip', win);
                % Is this the event video frame we've been waiting for?
                if (onsettime==-1 && pts >= timeOfEvent)
                    % Yes: This is the first frame with a pts timestamp that is
                    % equal or greater than the timeOfEvent, so 'vbl' is
                    % the exact time when the event was presented to the
                    % subject. Define it as onsettime:
                    onsettime = vbl;
                    
                    % Compare current pts to last one to see if the movie
                    % decoder skipped a frame at this crucial point in
                    % time. That would invalidate this trial.
                    if (pts - lastpts > 1.5*(1/fps))
                        % Difference to last frame is more than 1.5 times
                        % the expected difference under assumption 'no
                        % skip'. We skipped in the wrong moment!
                        rejecttrial=1;
                    end;
                end;
                
                % Keep track of the frames pts in order to check for skipped frames:
                lastpts=pts;
                
                % Delete the texture. We don't need it anymore:
                Screen('Close', movietexture);
                movietexture=0;
            end;
            
            % Done with drawing. Check the keyboard for subjects response:
            [keyIsDown, secs, keyCode]=KbCheck;
            if (keyIsDown==1)
                % Abort requested?
                if keyCode(esc)
                    % This signals abortion:
                    rejecttrial=-1;
                    % Break out of display loop:
                    break;
                end;
                
                % Space key pressed to indicate detection of event?
                if keyCode(space)
                    % Response too early (before event happened?)
                    if (onsettime==-1)
                        % Reject this trial:
                        rejecttrial=2;
                    else
                        % Valid response: Difference between 'secs' and
                        % 'onsettime' is the reaction time:
                        reactiontime=secs - onsettime;
                    end;
                end;
                
            end;
        end; % ...of display loop...
        
        % Stop movie playback, in case it isn't already stopped. We do this
        % by selection of a playback rate of zero: This will also return
        % the number of frames that had to be dropped to keep audio, video
        % and realtime in sync.
        droppedcount = Screen('PlayMovie', movie, 0, 0, 0);
        if (droppedcount > 0.2*framecount)
            % Over 20% of all frames skipped?!? Playback problems! We
            % reject this trial...
            rejecttrial=4;
        end;
        
        % Close the moviefile.
        Screen('CloseMovie', movie);
        
        % Check if aborted.
        if (rejecttrial==-1)
            % Break out of trial loop
            break;
        end;
        
        if (reactiontime==-1 && rejecttrial==0)
            rejecttrial=3;
        end;
        
        % Print out trials result if it was a valid trial:
        if (rejecttrial==0)
            fprintf('Trial %i valid: Reaction time was %f msecs.\n', i, 1000 * reactiontime);
        end;
        
        if (rejecttrial==1)
            fprintf('Trial %i rejected due to skip in video playback at time of event.\n', i);
        end;
        
        if (rejecttrial==2)
            fprintf('Trial %i rejected. False detection by subject.\n', i);
        end;
        
        if (rejecttrial==3)
            fprintf('Trial %i rejected. No detection by subject. Asleep?!?\n', i);
        end;
        
        if (rejecttrial==4)
            fprintf('Trial %i rejected. Way too many skips in movie playback!!!\n', i);
        end;
        
        % Wait for subject to release keys:
        KbReleaseWait;
        
    end; % Trial done. Next trial...
    
    % Done with the experiment. Close onscreen window and finish.
    ShowCursor;
    Screen('CloseAll');
    fprintf('Done. Bye!\n');
    return;
catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    sca;
    psychrethrow(psychlasterror);
end;
