function NewFastFilteredNoiseDemo(filtertype, rectSize, kwidth, scale, syncToVBL, dontclear, validate)
% NewFastFilteredNoiseDemo([filtertype=1][, rectSize=128][, kwidth=5][, scale=1][, syncToVBL=1][, dontclear=0][, validate=0])
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
    Screen('Preference', 'Enable3DGraphics', 0);
    Screen('Preference', 'Verbosity', 4);
    % Open fullscreen onscreen window on that screen. Background color is
    % gray, double buffering is enabled. Return a 'win'dowhandle and a
    % rectangle 'winRect' which defines the size of the window:
    [win, winRect] = Screen('OpenWindow', screenid, 128, [], [], [], [], [], mor(kPsychNeedFastBackingStore));
    
    % Build a filter kernel:
    stddev = kwidth / 3;

    switch(filtertype)
        case 0
            kernel = [];
        case 1
            kernel = fspecial('gaussian', kwidth, stddev);
        case 2
            kernel = fspecial('prewitt');
    end

    stype = 0;
    channels = 1;

    % Create processing proxy:
    filter = Screen('OpenProxy', win);

    % Enable it:
    Screen('HookFunction', filter, 'Enable', 'UserDefinedBlit');

    if filtertype > 0
        % Build shader from kernel:
        shader = EXPCreateStatic2DConvolutionShader(kernel, channels, 4, stype,1);
        
        % Attach our shader to proxies user defined processing chain:
        Screen('HookFunction', filter, 'PrependShader', 'UserDefinedBlit', '2DConvolutionOperator-1', shader);        
    else
        % Attach identity blit to proxies user defined processing chain:
        Screen('HookFunction', filter, 'PrependBuiltinFunction', 'UserDefinedBlit', 'Builtin:IdentityBlit', 'Builtin:IdentityBlit');
    end

    glFinish;
    
    % Compute destination rectangle locations for the random noise patches:
    
    % 'objRect' is a rectangle of the size 'rectSize' by 'rectSize' pixels of
    % our Matlab noise image matrix:
    objRect = SetRect(0,0, rectSize, rectSize);

    % ArrangeRects creates 'numRects' copies of 'objRect', all nicely
    % arranged / distributed in our window of size 'winRect':
    dstRect = ArrangeRects(numRects, objRect, winRect);

    % Now we rescale all rects: They are scaled in size by a factor 'scale':
    for i=1:numRects
        % Compute center position [xc,yc] of the i'th rectangle:
        [xc, yc] = RectCenter(dstRect(i,:));
        % Create a new rectange, centered at the same position, but 'scale'
        % times the size of our pixel noise matrix 'objRect':
        dstRect(i,:)=CenterRectOnPoint(objRect * scale, xc, yc);
    end

    %tex = Screen('OpenOffscreenWindow', win, 0, [0 0 rectSize rectSize]);            

    % Init framecounter to zero and take initial timestamp:
    count = 0; 
    glFinish;
    tstart = GetSecs;
    endtime = tstart + 5;
    ftex = 0;
    
    % Run noise image drawing loop for 20 seconds.
    while GetSecs < endtime
        % Generate and draw 'numRects' noise images:
        for i=1:numRects
            % Compute noiseimg noise image matrix with Matlab:
            % Normally distributed noise with mean 128 and stddev. 50, each
            % pixel computed independently:
            noiseimg=(50*randn(rectSize, rectSize) + 128);

            %noiseimg=zeros(rectSize, rectSize);
            %noiseimg(1:64,1)=255;
            

            if validate
                noiseimg=uint8(noiseimg);
            end
            
            % Convert it to a texture 'tex':
            tex=Screen('MakeTexture', win, noiseimg, [], [],[], 1);

            %tex=Screen('MakeTexture', win, transpose(flipud(noiseimg)), [], [],[], 2);
            %Screen('PutImage', tex, noiseimg);
            %Screen('DrawTexture', tex, intex, [], [], [], 0);
            
            % Draw the texture into the screen location defined by the
            % destination rectangle 'dstRect(i,:)'. If dstRect is bigger
            % than our noise image 'noiseimg', PTB will automatically
            % up-scale the noise image. We set the 'filterMode' flag for
            % drawing of the noise image to zero: This way the bilinear
            % filter gets disabled and replaced by standard nearest
            % neighbour filtering. This is important to preserve the
            % statistical independence of the noise pixels in the noise
            % texture! The default bilinear filtering would introduce local
            % correlatio
            if validate
                glFinish;
                tic
                % Transform texture into a new filtered texture:
                ftex = Screen('TransformTexture', tex, filter, [], ftex);
                
                % Draw it:
                Screen('DrawTexture', win, ftex, [], dstRect(i,:), [], 0);
                glFinish;
                gput = toc;
            else                
                % Transform texture into a new filtered texture:
                ftex = Screen('TransformTexture', tex, filter, [], ftex);
                
                % Draw it:
                Screen('DrawTexture', win, ftex, [], dstRect(i,:), [], 0);
            end
            
            if validate
                % Compute same convolution on CPU:
                noiseimg = single(noiseimg);
                tic
                ref = conv2(noiseimg, kernel, 'same');
                cput = toc;
                %ref = uint8(0.5 + ref);
            end
            
            % After drawing, we can discard the noise texture.
            Screen('Close', tex);
        end
        
        % Done with drawing the noise patches to the backbuffer: Initiate
        % buffer-swap. If 'asyncflag' is zero, buffer swap will be
        % synchronized to vertical retrace. If 'asyncflag' is 2, bufferswap
        % will happen immediately -- Only useful for benchmarking!
        Screen('Flip', win, 0, dontclear, asyncflag);
        
        if validate
            gpu = (Screen('GetImage', win, dstRect(1,:)));
            difference = double(gpu(:,:,1)) - double(ref);
            difference = difference(length(kernel):end-length(kernel), length(kernel):end-length(kernel));
            maxdiff = max(max(abs(difference)));
            fprintf('CPU %f secs, GPU %f secs, Speedup = %f, maxdiff = %f\n', cput, gput, cput / gput, maxdiff);
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

    if validate
        imagesc(difference);
        figure;
        imagesc(ref);
        figure;
        imagesc(gpu);
    end

    % Done. Close Screen, release all ressouces:
    sca;
catch
    % Our usual error handler: Close screen and then...
    sca;
    % ... rethrow the error.
    psychrethrow(psychlasterror);
end
