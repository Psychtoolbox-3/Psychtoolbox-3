function FastSeparableFilteredNoiseDemo(filtertype, rectSize, kwidth, scale, syncToVBL, dontclear, validate)
% FastSeparableFilteredNoiseDemo([filtertype=1][, rectSize=128][, kwidth=5][, scale=1][, syncToVBL=1][, dontclear=0][, validate=0])
%
% Demonstrates how to generate, filter and draw noise patches on-the-fly 
% in a fast way by use of GLSL fragment shaders.
% Use it to benchmark your system by varying the load. If you like this demo
% then also have a look at FastMaskedNoiseDemo that shows how to
% efficiently draw a masked stimulus by use of alpha-blending.
%
% filtertype = Type of filter to apply, see switch statement below for
% supported filters. Zero selects no filtering, 1 is Gaussian blur.
%
% rectSize = Size of the generated random noise image: rectSize by rectSize
%            pixels. This is also the size of the Psychtoolbox noise
%            texture.
%
% kwidth = For filters which support a varying kernel size, the kernel
% size. Will create a kwidth x kwidth convolution kernel.
%
% scale = Scalefactor to apply to texture during drawing: E.g. if you'd set
% scale = 2, then each noise pixel would be replicated to draw an image
% that is twice the width and height of the input noise image. In this
% demo, a nearest neighbour filter is applied, i.e., pixels are just
% replicated, not bilinearly filtered -- Important to preserve statistical
% independence of the random pixel values!
%
% syncToVBL = 1=Synchronize bufferswaps to retrace. 0=Swap immediately when
% drawing is finished. Value zero is useful for benchmarking the whole
% system, because your measured framerate will not be limited by the
% monitor refresh rate -- Gives you a feeling of how much headroom is left
% in your loop.
%
% dontclear = If set to 1 then the backbuffer is not automatically cleared
% to background color after a flip. Can save up to 1 millisecond on old
% graphics hardware.
%
% Example results on a Intel Pentium-4 3.2 Ghz machine with a NVidia
% GeForce 7800 GTX graphics card, running under M$-Windows XP SP3:
%
% Two patches, 256 by 256 noise pixels each, scaled by any factor between 1
% and 5 yields a redraw rate of 100 Hz.
%
% One patch, 256 by 256 noise pixels, scaled by any factor between 1
% and 5 yields a redraw rate of 196 Hz.
%
% Two patches, 128 by 128 noise pixels each, scaled by any factor between 1
% and 5 yields a redraw rate of 360 - 380 Hz.
% 
% One patch, 128 by 128 noise pixels, scaled by any factor between 1
% and 5 yields a redraw rate of 670 Hz.

% Abort script if it isn't executed on Psychtoolbox-3:
AssertOpenGL;

% Assign default values for all unspecified input parameters:
numRects = 1;

if nargin < 1 || isempty(filtertype)
    filtertype = 1; % Gaussian blur by default.
end

if nargin < 2 || isempty(rectSize)
    rectSize = 128; % Default patch size is 128 by 128 noisels.
end

if nargin < 3 || isempty(kwidth)
    kwidth = 5; % Kernel width is 5 x 5 by default.
end

if nargin < 4 || isempty(scale)
    scale = 1; % Don't up- or downscale patch by default.
end

if nargin < 5 || isempty(syncToVBL)
    syncToVBL = 1; % Synchronize to vertical retrace by default.
end

if syncToVBL > 0
    asyncflag = 0;
else
    asyncflag = 2;
end

if nargin < 6 || isempty(dontclear)
    dontclear = 0; % Clear backbuffer to background color by default after each bufferswap.
end

if dontclear > 0
    % A value of 2 will prevent any change to the backbuffer after a
    % bufferswap. In that case it is your responsibility to take care of
    % that, but you'll might save up to 1 millisecond.
    dontclear = 2;
end

if nargin < 7 || isempty(validate)
    validate = 0;
end

try
    % Find screen with maximal index:
    screenid = max(Screen('Screens'));

    % Initialize OpenGL. We need it for the image processing:
    InitializeMatlabOpenGL;
    
    % Open fullscreen onscreen window on that screen. Background color is
    % gray, double buffering is enabled. Return a 'win'dowhandle and a
    % rectangle 'winRect' which defines the size of the window:
    [win, winRect] = Screen('OpenWindow', screenid, 128);
    
    % Build a separable filter kernel:
    stddev = kwidth / 3;
stddev=30;

    switch(filtertype)
        case 0
            kernel = [];
        case 1
            % Standard 2D gauss kernel:
            kernel  = fspecial('gaussian', kwidth, stddev);
            % Same kernel split up into two 1D kernels:
            kernel1 = fspecial('gaussian', [kwidth, 1], stddev);
            kernel2 = fspecial('gaussian', [1, kwidth], stddev);            
        case 2
            kernel = fspecial('prewitt');
    end
    
    stype = 2;
    channels = 1;

    if filtertype > 0
        % Build shaders from kernels:
        shader1 = EXPCreateStatic2DConvolutionShader(kernel1, channels, 1, stype, 1);
        shader2 = EXPCreateStatic2DConvolutionShader(kernel2, channels, 1, stype, 1);
    end

    % Create a floating-point precision FBO for intermediate result of first filter pass:
    [convolvefbo, convolvetex] = moglCreateFBO(rectSize, rectSize, 1, 4);
    
    % Sync us to the GPU so our timing benchmarks are ok:
    glFinish;
    
    % Compute destination rectangle locations for the random noise patches:    
    % 'objRect' is a rectangle of the size 'rectSize' by 'rectSize' pixels of
    % our Matlab noise image matrix:
    objRect = SetRect(0,0, rectSize, rectSize);
    % dstRect = CenterRect(objRect, winRect);
    dstRect = objRect;

    % Init framecounter to zero and take initial timestamp:
    count = 0; 
    glFinish;
    tstart = GetSecs;
    endtime = tstart + 15;
    
    % Run noise image drawing loop for 20 seconds.
    while GetSecs < endtime
        % Generate and draw 'numRects' noise images:
        for i=1:numRects
            % Compute noiseimg noise image matrix with Matlab:
            % Normally distributed noise with mean 128 and stddev. 50, each
            % pixel computed independently:
            if (count==0 || validate < 2) noiseimg=(50*randn(rectSize, rectSize) + 128); end

            % Validation enabled?
            if validate == 1
                % Cast it to uint8 so GPU and CPU get exactly the same
                % input data:
                noiseimg=uint8(noiseimg);
            end
            
            % Convert it to a texture 'tex':
            tex=Screen('MakeTexture', win, noiseimg);
            intex=Screen('GetOpenGLTexture', win, tex);
            cotex=Screen('SetOpenGLTexture', win, [], convolvetex, GL.TEXTURE_RECTANGLE_EXT);
            
            if validate == 1
                glFinish;
                tic
            end

            % Pass 1: Blit texture into our FBO, apply first kernel:
            %moglChooseFBO(convolvefbo);
            moglChooseFBO(0);
            glUseProgram(shader1);
            %moglBlitTexture(intex, [], [], [], [], 0, GL.CLAMP_TO_BORDER);
            Screen('DrawTexture', win, tex, [], dstRect, [], 0);
            % Pass 2: Blit intermediate result into framebuffer, apply
            % second kernel:
            %moglChooseFBO(0);
            %glUseProgram(shader2);
            %Screen('DrawTexture', win, cotex, [], dstRect, [], 0);
            %moglBlitTexture(convolvetex, [], [], [], [], 0, GL.CLAMP_TO_BORDER);

            % Done. Disable shader:
            glUseProgram(0);

            if validate == 1
                glFinish;
                gpu = toc
            end
            
            if validate == 1
                % Compute same convolution on CPU:
                noiseimg = single(noiseimg);
                tic
                ref = conv2(noiseimg, single(kernel1), 'same');
                %ref = conv2(single(ref), single(kernel2), 'same');
                %ref = conv2(noiseimg, (kernel), 'same');
                cpu = toc
                ref = uint8(0.5 + ref);
            end
            
            % After drawing, we can discard the noise texture.
            Screen('Close', tex);
        end
        
        % Done with drawing the noise patches to the backbuffer: Initiate
        % buffer-swap. If 'asyncflag' is zero, buffer swap will be
        % synchronized to vertical retrace. If 'asyncflag' is 2, bufferswap
        % will happen immediately -- Only useful for benchmarking!
        Screen('Flip', win, 0, dontclear, asyncflag);
        
        if validate == 1
            gpu = Screen('GetImage', win, dstRect);
            difference = gpu(:,:,1) - ref;
            difference = difference(length(kernel):end-length(kernel), length(kernel):end-length(kernel));
            maxdiff = max(max(difference))
        end
        
        % Increase our frame counter:
        count = count + 1;
    end

    % We're done: Output average framerate:
    glFinish;
    telapsed = GetSecs - tstart
    updaterate = count / telapsed
    
    % Disable shader: Standard fixed-function pipeline is activated.
    glUseProgram(0);

    % Delete our FBO and its associated color buffer texture:
    moglDeleteFBO(convolvefbo);
    
    if validate == 1
        imagesc(difference);
        figure;
        imagesc(noiseimg);
    end

    % Done. Close Screen, release all ressouces:
    Screen('CloseAll');
catch
    % Our usual error handler: Close screen and then...
    Screen('CloseAll');
    % ... rethrow the error.
    psychrethrow(psychlasterror);
end
