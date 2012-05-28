function DrawingSpeedTest(n,primitivetype,mode,gpumeasure)
% DrawingSpeedTest([n=800][,primitivetype=0][,mode=0][,gpumeasure=0])
%
% Tests batch-drawing performance of some Screen functions. Batch drawing
% is a way to submit multiple primitives, e.g., Filled Rects, at once. This
% way the Matlab call overhead and some internal setup overhead is saved
% and drawing should be significantly faster than when issuing single
% drawing commands in a loop.
%
% This currently only tests filled rects and framed rects as well as filled
% ovals. It also provides a way to test drawing by texture mapping.
% Dots and Lines are nicely demonstrated by DotDemo and LinesDemo.
%
% The optional parameter n allows to specifiy the number of primitives to
% draw each frame, default is 800. The test loop will draw 1000 identical
% frames and measure the time needed.
%
% 'primitivetype' type of primitive: 0 = filled rects, 1 = framed rects, 2
% = filled ovals, 3 = framed ovals.
%
% 'mode' type of drawing: 0 = One by one submission (slowest), 1 = batch
% submission, 2 = texture mapping for drawing, 3 = texture mapping with
% batch submission of textures.
%
% 'gpumeasure' measurement type: 0 = Only on cpu. 1 = Measure exact
% execution time of drawing commands on GPU's that support this feature.
% Plot result at end of run.
%

% History:
% 04/30/07 mk Wrote it.
% 12/21/09 mk Add support for GPU measurement via EXT_timer_query on
%             supported hardware.

if nargin < 1 || isempty(n)
    n = 800;
end

if nargin < 2 || isempty(primitivetype)
    primitivetype = 0;
end

if nargin < 3 || isempty(mode)
    mode = 0;
end

if nargin < 4 || isempty(gpumeasure)
    gpumeasure = 0;
end

% Check proper PTB installation:
AssertOpenGL;

% Open window with gray background on secondary display (if any):
screenid = max(Screen('Screens'));
[win winrect] = Screen('OpenWindow', screenid, 128, [], [], [], [], 0);

% Setup stim parameters:
w=RectWidth(winrect);
h=RectHeight(winrect);
sizeX=80;
sizeY=80;
msize = min(w,h);

% In mode 2 or 3, we draw a prototype stimulus into an offscreen window,
% using the alpha-channel to mask out all non-stimulus shape pixels. We
% will use this texture image of a prototype stim later on to draw the
% stimulus shape, rescaling it and changing its color and alpha as
% requested:
if mode == 2 || mode == 3
    % Need alpha-blending - Enable it:
    Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    % Open offscreen window 'template' of reasonable size msize x msize
    % pixels. Set it to a black transparent background [r g b a] = [0 0 0 0]
    template = Screen('OpenOffscreenWindow', win, [0 0 0 0], [0 0 msize msize]);
    % Draw prototype primitive which covers the full window and has a color
    % and alpha value of 255 ie. white and fully opaque:
    switch primitivetype
        case 0,
            Screen('FillRect', template, 255);
        case 1,
            Screen('FrameRect', template, 255, [0 0 msize msize], 50);
        case 2,
            Screen('FillOval', template, 255);
        case 3,
            Screen('FrameOval', template, 255);
    end
end

% Generate a matrix which specs n filled rectangles, with randomized
% position, color and (dot-,line-)size parameter
for j = 1:n
    posx = rand * w;
    posy = rand * h;
    colors(j, 1:3) = [ rand * 255, rand * 255, rand * 255];
    myrect(j, 1:4) = floor([ posx, posy, posx + rand * sizeX, posy + rand * sizeY]);
    sizes(j) = floor(rand * 10)+1;
end

% Need color and shape vectors transposed, i.e., 4 row by n columns
% matrices:
colors = transpose(colors);
myrect = transpose(myrect);

gpudur = zeros(1, 1000);

% Initially sync us to retrace, take start time t1:
t1 = Screen('Flip', win);

% Timing loop, 1000 trials:
for i=1:1000
    if gpumeasure
        % Start GPU timer: gpumeasure will be true if this
        % is actually supported and will return valid results:
        gpumeasure = Screen('GetWindowInfo', win, 5);
    end

    % Batch draw:
    if mode < 2
        switch primitivetype
            case 0,
                if mode == 1
                    % Batch drawing version of FillRect - submit all
                    % primitives at once:
                    Screen('FillRect', win, colors, myrect);
                end

                if mode == 0
                    for j=1:n
                        % One part at a time submission:
                        Screen('FillRect', win, colors(:,j)', myrect(:,j)');
                    end
                end
            case 1,
                if mode == 1
                    Screen('FrameRect', win, colors, myrect, sizes);
                end

                if mode == 0
                    for j=1:n
                        Screen('FrameRect', win, colors(:,j)', myrect(:,j)', sizes(j));
                    end
                end

            case 2,
                if mode == 1
                    Screen('FillOval', win, colors, myrect);
                end

                if mode == 0
                    for j=1:n
                        Screen('FillOval', win, colors(:,j)', myrect(:,j)');
                    end
                end
            case 3,
                if mode == 1
                    Screen('FrameOval', win, colors, myrect, sizes);
                end

                if mode == 0
                    for j=1:n
                        Screen('FrameOval', win, colors(:,j)', myrect(:,j)', sizes(j));
                    end
                end

        end
    else
        % mode > 2: Use DrawTexture to draw template primitive:
        if mode == 2
            for j=1:n
                % DrawTexture in a loop:
                Screen('DrawTexture', win, template, [], myrect(:,j)', [], 0, [], colors(:,j)');
            end;
        else
            % Batch drawing version DrawTextures:
            Screen('DrawTextures', win, template, [], myrect, 0, 0, 0, colors);
        end
    end

    % Flip it. Don't clear buffers, don't sync to retrace. We want the raw
    % speed, nothing limited by monitor refresh:
    Screen('Flip', win, 0, 2, 2);

    % Result of GPU time measurement expected?
    if gpumeasure
        % Retrieve results from GPU load measurement:
        % Need to poll, as this is asynchronous and non-blocking,
        % so may return a zero time value at first invocation(s),
        % depending on how deep the rendering pipeline is:
        while 1
            winfo = Screen('GetWindowInfo', win);
            if winfo.GPULastFrameRenderTime > 0
                break;
            end
        end

        % Store it:
        gpudur(i) = winfo.GPULastFrameRenderTime;
    end

    % Next iteration.
end

% Make sure the GPU is idle:
Screen('DrawingFinished', win, 2, 1);

% Take end time, compute and print the stats:
telapsed = GetSecs - t1;
tavg = telapsed;
tperrect = tavg / n;

fprintf('Rendered 1000 frames, each with %i primitives of size %i x %i.\n', n, sizeX,sizeY);
fprintf('Total cpu time %6.6f seconds. Time per primitive %6.6f msecs.\n', telapsed, tperrect);

%Done.
Screen('CloseAll');

if any(gpudur)
    gpudur = 1000 * gpudur;
    plot(gpudur);
    title('Drawtime in msecs per frame');
    fprintf('Mean drawtime on GPU is %f msecs per frame, stddev = %f msecs, median %f msecs, time per primitive %6.6f msecs.\n', mean(gpudur), std(gpudur), median(gpudur), mean(gpudur)/n);
end

return;
