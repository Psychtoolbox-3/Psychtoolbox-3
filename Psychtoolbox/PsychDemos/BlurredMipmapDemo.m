function BlurredMipmapDemo(deviceIndexOrMoviename, gazeRadius, customMipmap)
% BlurredMipmapDemo - Show blurring of live video or movies via adaptive mipmapping.
%
% This demo shows how to selectively blur regions of a live video stream
% captured from an attached camera, or a live playing movie. The method
% employed here uses a GLSL shader to select an individual blur level for
% each output pixel. Pixels close to the mouse cursor (which simulates gaze
% position, e.g., as aquired by an eyetracker) will be drawn sharp, whereas
% pixels at larger radial distance will be drawn increasingly unsharp by
% low-pass filtering them.
%
% The method of blurring here is computationally efficient and fast,
% because it generates and uses a image resolution pyramid to precompute
% different versions of the video image at different resolutions. Then the
% shader just looks up color information for each individual output pixel
% from the proper "level of detail" in the pyramid of low pass filtered
% images. The resolution pyramid is auto-generated in a fast, GPU
% accelerated way for each image by use of OpenGL "mip-mapping". This
% method uses a fixed low-pass filter for downsampling, usually a box
% filter, but the exact type of filter used is system dependent and can
% vary across different graphics cards, graphics driver versions or
% operating systems.
%
% The demo requires recent graphics hardware and won't work on old graphics
% cards.
%
% Usage:
%
% BlurredMipmapDemo([deviceIndexOrMoviename=0][, gazeRadius=25][, customMipmap=0])
%
% 'gazeRadius' Some falloff radius in pixels. Controls how quickly the
% resolution degrades with increasing distance from the cursorposition.
% Smaller values give a smaller simulated "foveal area".
%
% 'deviceIndexOrMoviename' = Optional: Either the index of the video
% capture device, or the full path to a movie file. Defaults to
% auto-selected default video source.
%
% 'customMipmap' Optional: If non-zero, use a custom shader for
% downsampling during mipmap building, instead of the relatively simple
% builtin filter of the GPU.
%
% Press any key to finish the demo.
%

% History:
% 28.08.2012  mk  Written.

% Child protection:
AssertOpenGL;

% Select maximum display id for output:
screenId = max(Screen('Screens'));

if nargin < 1 || isempty(deviceIndexOrMoviename)
    deviceIndexOrMoviename = [];
end

if nargin < 2 || isempty(gazeRadius)
    gazeRadius = 25;
end

if nargin < 3 || isempty(customMipmap)
    customMipmap = 0;
end

% Is it a movie file name?
if ischar(deviceIndexOrMoviename)
    % Yes: Playback the movie.
    moviename = deviceIndexOrMoviename;
    ismovie = 1;
else
    % No: Nothing or a numeric video capture deviceindex. Do videocapture.
    ismovie = 0;
end

if customMipmap == 0
    % Use automatic mipmap generation by GPU driver builtin method:
    dontautomip = 0;
else
    % Don't generate mipmaps automatically, roll our own downsampling:
    if ismovie
        % This is the correct flag for 'GetMovieImage':
        dontautomip = 8; 
    else
        % This is the correct flag for 'GetCapturedImage':
        dontautomip = 16;
    end
end

try
    % Open onscreen window with black background color:
    win = Screen('OpenWindow', screenId, 0);
    [w, h] = Screen('WindowSize', win); %#ok<*ASGLU>
    
    % Load & Create a GLSL shader for adaptive mipmap lookup:
    shaderpath = fileparts(mfilename('fullpath'));
    shader = LoadGLSLProgramFromFiles([shaderpath filesep 'BlurredMipmapDemoShader'], 1);
    
    % Bind texture unit 0 to shader as input source for the mip-mapped video image:
    glUseProgram(shader);
    glUniform1i(glGetUniformLocation(shader, 'Image'), 0);
    glUseProgram(0);

    % Load & Create a GLSL shader for downsampling during MipMap pyramid
    % creation. This specific example shader uses a 3-by-3 gaussian filter
    % kernel with a standard deviation of 1.0:
    mipmapshader = LoadGLSLProgramFromFiles([shaderpath filesep 'MipMapDownsamplingShader.frag.txt'], 1);
    
    if ismovie
        % Open movie and start its playback:
        movie = Screen('OpenMovie', win, moviename);
        Screen('PlayMovie', movie, 1, 1, 1);
    else
        % Open video capture device 'deviceIndex':
        grabber = Screen('OpenVideoCapture', win, deviceIndexOrMoviename);
        
        % Start video capture:
        Screen('StartVideoCapture', grabber, realmax, 1);
    end
    
    % Change cursor to a cross-hair:
    ShowCursor('CrossHair');
    
    tex = [];
    count = 0;
    t=GetSecs;
    
    % Repeat until keypress or timeout of 10 minutes:
    while ((GetSecs - t) < 600) && ~KbCheck
        % Wait for next video image or movie frame. The first '1' asks to
        % wait blocking for the next image, instead of polling. 'tex' is
        % returned for recycling to increase efficiency. The specialMode
        % flag '1' asks Screen() to return the texture as a GL_TEXTURE_2D
        % texture. This is needed for the automatic mip-mapping for
        % low-pass filtering to work:
        if ismovie
            % Get next movie frame:
            tex = Screen('GetMovieImage', win, movie, 1, [], 1 + dontautomip);
        else
            % Get next video frame:
            tex = Screen('GetCapturedImage', win, grabber, 1, tex, 1 + dontautomip);
        end
        
        % Valid 'tex'ture with new videoframe returned?
        if tex > 0
            % Yes. Get current "center of gaze" as simulated by the current
            % mouse cursor position:
            [gazeX, gazeY] = GetMouse(win);
            
            % Flip y-axis direction -- Shader has origin bottom-left, not
            % top-left as GetMouse():
            gazeY = h - gazeY;
            
            % Need to create mipmaps ourselves?
            if dontautomip
                % Yes: Use PTB's utility function for shader-based mipmap
                % generation with more control, but less performance. We
                % use 'mipmapshader' as downsampling shader, and '1' enable
                % bilinear filtering for the shader:
                CreateResolutionPyramid(tex, mipmapshader, 1);
            end
            
            % Draw new video texture from framegrabber. Apply GLSL 'shader'
            % during drawing. Use filterMode 3: This will automatically
            % generate the mipmap image resolution pyramid, then use the
            % 'shader' to adaptively lookup lowpass filtered pixels from
            % the different blur-levels of the image pyramid to simulate
            % typical foveation effect of decreasing resolution with
            % increasing distance to center of fixation. We pass gaze
            % center (gazeX, gazeY) and radius of foveation 'gazeRadius' as
            % auxParameters. auxParameters must always have a multiple of 4
            % components, so we add a zero value for padding to length 4:
            Screen('DrawTexture', win, tex, [], [], [], 3, [], [], shader, [], [gazeX, gazeY, gazeRadius, 0]);
            
            % Need to close texture if movie is played back:
            if ismovie
                Screen('Close', tex);
            end
            
            % Show it at next video refresh:
            Screen('Flip', win);
        end
        
        % Increase frame counter:
        count = count + 1;
    end
    
    % Display some timing stats:
    telapsed = GetSecs - t;
    avgfps = count / telapsed;
    fprintf('Average frames per second: %f.\n', avgfps);
    
    if ismovie
        % Stop and close movie:
        Screen('PlayMovie', movie, 0);
        Screen('CloseMovie', movie);
    else
        % Stop and close video source:
        Screen('StopVideoCapture', grabber);
        Screen('CloseVideoCapture', grabber);
    end
    
    % Close down everything else:
    Screen('CloseAll');
    
catch %#ok<*CTCH>
    % Error handling, emergency shutdown:
    sca;
    psychrethrow(psychlasterror);
end
