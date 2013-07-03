function LoadMovieIntoTexturesDemo(moviename, fromTime, toTime, indexisFrames, benchmark, async, preloadSecs, specialflags, pixelFormat)
%
% LoadMovieIntoTexturesDemo(moviename [, fromTime=0][, toTime=end][, indexisFrames=0][, benchmark=0][, async=0][, preloadSecs=1][, specialflags=0][, pixelFormat=4])
%
% A demo implementation on how to load a movie into standard
% Psychtoolbox textures for precisely controlled presentation timing and
% presentation order.
%
% Parameters:
% moviename - Filename of moviefile to use. If none is provided, then the
% simple DualDiscs collision movie is used which is part of PTB.
%
% fromTime - Start time (in seconds) from which the movie should be read
% into textures. Defaults to start of movie, if not provided.
%
% toTime - End time (in seconds) upto which the movie should be read
% into textures. Defaults to end of movie, if not provided.
%
% indexIsFrames - If set to 1 then the fromTime and toTime parameters are
% interpreted as frameindex (starting with 0 for the first frame), instead
% of seconds. If set to 2 then presentation timestamps are ignored for the
% decision when to stop loading a movie. This is needed for certain movies
% with a broken encoding of time, as those would stop loading way too early
% otherwise.
%
% benchmark - If you set this parameter to 1, the demo will compute the
% time it takes to load the movie, and - after the movie has been loaded
% into textures - the maximum display speed without syncing to vertical
% refresh. All visual progress feedback is disabled. This is mostly
% useful to benchmark different movie file formats and codecs for their
% relative efficiency on a given machine. A setting of 2 will additionally
% disable test of abort keys during loading of the movie - to discount
% possible pollution of the benchmark results with time spent checking the
% keyboard. A setting of 3 will also skip keyboard queries during the load
% phase, additionally it will not queue up textures in video memory, but
% instead load a texture and then immediately discard it. This is more
% resembling typical movie playback, where a frame is pulled, presented on
% screen, then discarded. It tests a different decoding path in Screen(). A
% setting of 4 will only test decoding performance of the movie, but omit
% creation of actual Psychtoolbox textures for presentation. This allows to
% separate computation time spent in the video decoder from time consumed
% by the graphics driver or graphics card.
%
% async - If you set this parameter to 4, the video decoding engine will
% prebuffer frames ahead of time, up to 'preloadSecs' seconds worth of
% video data. This is unsuitable for playback with audio-video sync, but
% for pure video playback it can decouple decoding from presentation
% further and provide a potential performance boost for very demanding
% playback scenarios.
%
% preloadSecs - How many seconds of video to prebuffer if async == 4?
% Specify a maximum amount in seconds, or the value -1 for unlimited
% prebuffering.
%
% specialflags - Special flags for 'OpenMovie'. E.g., a setting of 1 will
% try to use YUV textures for higher performance, if the installed graphics
% card supports this. A setting of 4 will always use YUV decoding, by use
% of out own builtin YUV decoder, which may be more limited in
% functionality and flexibility, but helpful if highest performance is a
% requirement.
%
% pixelFormat - Format of video texture to create: 1 = Luminance/Grayscale
% only, 2 = Luminance+Alpha, 3 = RGB color, 4 = RGBA, 5 = YUV-422 packed
% pixel, 6 = YUV-I420 planar format. 5 and 6 = Y8-Y800 planar luminance
% only format. Not all formats are supported by all GPU's, operating
% systems and video codecs. Defaults to 4 = RGBA 8 Bit per color channel.
%
%
% How the demo works: Read the source code - its well documented ;-)
%
% This demo "preloads" the movie into textures:
% The whole movie gets read into PTB textures before start of trial. Then
% you show the textures in quick succession like in PlayMoviesDemo. The
% advantage of this approach is exact control over display timing and
% display order: You can show frames in any order you like at any rate you
% like (and that your hardware likes). Disadvantage: Longer trial
% setup time for loading the whole movie, higher memory consumption (keep
% all n frames in memory instead of only the current one) and the inability
% to play sound in sync with video.
%
% To make it a bit more interesting, you can use this demo to "browse" short
% videoclips, e.g., for selection of interesting parts...
%

% History:
% 12/25/05  mk  Wrote it.
% 02/03/06  mk  Adapted for use on Windows.
% 09/03/09  mk  Add support for frameindex seeking instead of timeindex.
% 30/06/12  mk  Add benchmarking support for new GStreamer modes.
% 17/06/13  mk  Cleanup for QT removal.

% Child protection: Make sure we run on the OSX / OpenGL Psychtoolbox.
% Abort if we don't:
AssertOpenGL;

% Switch KbName into unified mode: It will use the names of the OS-X
% platform on all platforms in order to make this script portable:
KbName('UnifyKeyNames');

esc=KbName('ESCAPE');
space=KbName('SPACE');
right=KbName('RightArrow');
left=KbName('LeftArrow');

if nargin < 1
    moviename = [];
end

if isempty(moviename)
    % Default movie is our own disc collision movie:
    moviename = [ PsychtoolboxRoot 'PsychDemos/MovieDemos/DualDiscs.mov' ];
end;

if nargin < 2
    % Default to beginning of movie:
    fromTime=0;
end;

if isempty(fromTime)
    fromTime=0;
end

if nargin < 3
    % Default to end of movie:
    toTime=100000;
end;

if isempty(toTime)
    % Default to end of movie:
    toTime=100000;
end;

if nargin < 4
    indexisFrames = [];
end

if isempty(indexisFrames)
    indexisFrames = 0;
end

if nargin < 5
    benchmark = [];
end;

if isempty(benchmark)
    benchmark = 0;
end;

if nargin < 6
    async = [];
end

if nargin < 7
    preloadSecs = [];
end

if nargin < 8
    specialflags = [];
end

if nargin < 9
    pixelFormat = [];
end

fprintf('Loading movie %s ...\n', moviename);
try    
    % Background color will be a grey one:
    background=[128, 128, 128];

    % Open onscreen window. We use the display with the highest number on
    % multi-display setups:
    screen=max(Screen('Screens'));

    % This will open a screen with default settings, aka black background,
    % fullscreen, double buffered with 32 bits color depth:
    [win, winrect] = Screen('OpenWindow', screen); % , 0, [0 0 800 600]);
    
    % Hide the mouse cursor:
    HideCursor;
    
    % Clear screen to background color:
    Screen('FillRect', win, background);
    
    % Show instructions...
    tsize=20;
    Screen('TextSize', win, tsize);
    Screen('DrawText', win, 'Loading movie into textures. Please wait...',40, 100);    

    % Flip to show the startup screen:
    Screen('Flip',win);
        
    % Open the moviefile and query some infos like duration, framerate,
    % width and height of video frames...
    
    % Some legacy Screen() mex files don't support the pixelFormat flag,
    % therefore we try to do without it if it is not used anyway, avoiding
    % an error abort:
    if ~isempty(pixelFormat)
        % Use pixelFormat:
        [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename, async, preloadSecs, specialflags, pixelFormat);
    else
        % Legacy compatible call without pixelFormat:
        [movie movieduration fps imgw imgh] = Screen('OpenMovie', win, moviename, async, preloadSecs, specialflags);
    end
    
    % Move to requested timeindex where texture loading should start:
    if indexisFrames
        Screen('SetMovieTimeIndex', movie, fromTime, 1);
    else
        Screen('SetMovieTimeIndex', movie, fromTime, 0);
    end
    
    movietexture=0;     % Texture handle for the current movie frame.
    lastpts=-1;         % Presentation timestamp of last frame.
    pts=-1;
    count=0;            % Number of loaded movie frames.
    
    % Preallocate texids and texpts if maximum possible allocation size is
    % known:
    if (movieduration > 0) && (movieduration < intmax) && (fps > 0)
        texids = zeros(1, ceil(movieduration * fps));
        texpts = zeros(1, ceil(movieduration * fps));
    end
    
    % In benchmark mode 4 we don't actually convert decoded video frames
    % into textures. Instead we immediately discard them, so we just get a
    % benchmark for how long actual video decoding takes:
    if benchmark == 4
        % Discard fetched video frames before texture creation:
        dontfetch = 2;
    else
        % Convert decoded video frames into textures:
        dontfetch = 0;
    end

    % In async buffering mode we have to start the playback engine, so that
    % any ahead-of-time predecoding and prebuffering is actually happening:
    if bitand(async, 4)
        % We run playback at 100x the normal speed, non-looped, without
        % sound:
        Screen('PlayMovie', movie, 100, 0, 0);
        WaitSecs(5);
    end
    
    tloadstart=GetSecs;
    
    % Movie to texture conversion loop:
    while (movietexture>=0) && ((~indexisFrames && (pts < toTime)) || (indexisFrames && (fromTime + count <= toTime)))
            % This call waits for arrival of a new frame from the movie. If
            % a new frame is ready, it converts the video frame into a
            % Psychtoolbox texture image and returns a handle in
            % 'movietexture'. 'pts' contains a so called presentation
            % timestamp. That is the time (in seconds since start of movie)
            % at which this video frame should be shown on the screen.
            % In case that no valid video texture is available, this function
            % will return a zero texture handle to indicate this. If no new
            % texture will become available anymore, because the end of the
            % movie is reached, it will return a handle of -1 to indicate end
            % of movie.
            % The 1 - flag means: Wait for arrival of new frame. Next
            % invocation of this command will retrieve the next frame in
            % the movie.
            % A 'dontfetch' setting of 2 would just discard the new frame,
            % instead of returning a real PTB texture it returns a value of
            % 1 to indicate discard of the frame. This can be used to
            % quickly skip through movies by skipping frames, and for
            % benchmarking of the video decoder without the overhead of
            % actual texture creation and drawing:
            [movietexture pts] = Screen('GetMovieImage', win, movie, 1, [], [], dontfetch);
            
            % Valid and *new* texture? If pts would be *smaller* than
            % lastpts then we would have ran over the end of movie - in
            % that case, the time will automatically wrap around to zero.
            % If we don't check for this, we'll have an infinite loop!
            if (movietexture > 0) && ((pts >= lastpts) || (indexisFrames == 2))
                % Store its texture handle and exact movie timestamp in
                % arrays for later use:
                count=count + 1;
                
                % Only store textures for benchmark modes 0 to 3.
                if (benchmark ~= 3) && (benchmark ~= 4)
                    texids(count)=movietexture;
                    texpts(count)=pts;
                else
                    % Mode 3 or 4. Draw and delete texture immediately,
                    % resembling typical movie playback, or - in mode 4 -
                    % don't draw at all:
                    if benchmark ~=4
                        %Screen('TextSize', movietexture, 60);
                        %DrawFormattedText(movietexture, num2str(count), 'center', 'center', [255 255 0]);
                        Screen('DrawTexture', win, movietexture, [], winrect, [], 0);
                        Screen('Close', movietexture);
                        Screen('Flip', win, [], 2, 2);
                    end
                end
                lastpts=pts;
            else
                break;
            end

            % Allow for abortion...
            if (benchmark < 2) && KbCheck
                break;
            end
            
            if (benchmark==0)
                % Show the progress text:
                Screen('DrawText', win, ['Loaded texture ' num2str(count) '...'],40, 100);
                Screen('Flip',win);
            end
    end;
    
    if (benchmark>0)
        % Compute movie load & conversion rate in frames per second.
        loadrate = count / (GetSecs - tloadstart);
        % Compute same rate in Megapixels per second:
        loadvolume = loadrate * imgw * imgh / 1024 / 1024;

        fprintf('Movie to texture conversion speed is %f frames per second == %f Megapixels/second.\n', loadrate, loadvolume);
    end;

    if bitand(async, 4)
        % Stop "playback":
        Screen('PlayMovie', movie, 0, 0, 0);
    end
    
    if (benchmark == 3) || (benchmark == 4)
        count = 0;
    end
    
    % Ok, now the requested part of the movie has been (hopefully) loaded
    % and converted into standard PTB textures. We can simply use the
    % 'DrawTexture' command in a loop to show the textures...
    
    % Clear screen to background color:
    Screen('Flip',win);
    
    currentindex=1;
    autoplay=0;
    
    totalcount = 0;
    tstart=GetSecs;
    
    % Browse and Draw loop:
    while(count>0)
        % Draw texture 'currentindex'
        Screen('DrawTexture', win, texids(currentindex));

        if (benchmark==0)
            % Draw some help text:
            [x, y]=Screen('DrawText', win, 'Press left-/right cursor key to navigate in movie, SPACE to toggle playback, ESC to exit.',10, 40); %#ok<ASGLU>
            [x, y]=Screen('DrawText', win, ['Framerate(fps): ' num2str(fps) ', total duration of movie (secs) ' num2str(movieduration)],10, y+10+tsize); %#ok<ASGLU>

            % Draw info on current position in movie:
            Screen('DrawText', win, ['Frame ' num2str(currentindex) ' of ' num2str(count) ' : Timeindex(secs) = ' num2str(texpts(currentindex))], 10, y + 10 + tsize);

            % Show drawn stuff:
            Screen('Flip', win);

            % Check for key press:
            [keyIsDown, secs, keyCode]=KbCheck; %#ok<ASGLU>
            if keyIsDown
                if (keyCode(esc))
                    % Exit
                    break;
                end;

                if (keyCode(space))
                    % Toggle playback on space.
                    autoplay=1-autoplay;
                end;

                if (keyCode(right) && currentindex<count)
                    % One frame forward:
                    currentindex=currentindex+1;
                end;
                if (keyCode(left) && currentindex>1)
                    % One frame backward:
                    currentindex=currentindex-1;
                end;

                % Wait for key-release:
                KbReleaseWait;
            end;
        else
            % Benchmark mode: Repeat until Keypress.
            Screen('Flip', win, 0, 2, 2);
            totalcount=totalcount + 1;
            autoplay=1;
            
            if KbCheck
                break;
            end;
        end;
        
        % Update frameindex if in autoplay mode:
        if autoplay==1
            currentindex=mod(currentindex, count) + 1;
        end;
    end;
    
    % Done. Flip a last time to show grey background:
    Screen('Flip', win);
    
    if (benchmark>0)
        playbackrate = totalcount / (GetSecs - tstart);
        playbackvolume = playbackrate * imgw * imgh / 1024 / 1024;
        fprintf('Movietexture playback rate is %f frames per second == %f Megapixels/second.\n', playbackrate, playbackvolume);
    end;
    
    if benchmark ~= 3
        % This will release all textures...
        Screen('Close', texids(1:count));
    end
    
    % Close movie file.
    Screen('CloseMovie', movie);
   
    ShowCursor;
    Screen('CloseAll');
    fprintf('Done. Bye!\n');
    return;

catch %#ok<CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    ShowCursor;
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end;
