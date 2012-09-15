function BlurredVideoCaptureDemo(filtertype, kwidth, deviceIndex)
% BlurredVideoCaptureDemo - Show application of image filters to videocaptured images.
%
% This demo shows how to apply 2D convolution kernels to a live video
% stream captured from an attached default camera. Convolution is not
% performed on the cpu by Matlab/Octave, but executed hardware-accelerated
% on the GPU for potentially much higher speed.
%
% Usage:
%
% BlurredVideoCaptureDemo([filtertype=1][, kwidth=11][, deviceIndex=0])
%
% 'filtertype' is optional and can be:
%
% 0 = None.
% 1 = Gaussian blur - This is the default.
% 2 = Prewitt operator.
% 3 = Unsharp operator.
% 4 = Sobel operator.
% 5 = Log operator.
%
% The filters are selected via the fspecial() function, and converted into
% proper convolution kernels. You may need the Matlab image processing
% toolobox or equivalent Octave package installed, for fspecial() to be
% available.
%
% 'kwidth' = Optional convolution kernel width in pixels. Defaults to 11.
%
% 'deviceIndex' = Optional index of the video capture device. Defaults to
% auto-selected default video source.
%
% Press any key to finish the demo.
%

AssertOpenGL;
screen=max(Screen('Screens'));

if nargin < 1 || isempty(filtertype)
    filtertype = 1;
end

if nargin < 2 || isempty(kwidth)
    kwidth=11;
end

if nargin < 3 || isempty(deviceIndex)
    deviceIndex = [];
end

try
    InitializeMatlabOpenGL;
    
    win=Screen('OpenWindow', screen, 0);
    
    % Initial flip to a blank screen:
    Screen('Flip',win);
    
    % Set text size for info text. 24 pixels is also good for Linux.
    Screen('TextSize', win, 24);
    
    % Build a filter kernel:
    stddev = kwidth / 3;
    
    switch(filtertype)
        case 0
            kernel = [];
        case 1
            kernel = fspecial('gaussian', kwidth, stddev);
        case 2
            kernel = fspecial('prewitt');
        case 3
            kernel = fspecial('unsharp');
        case 4
            kernel = fspecial('sobel');
        case 5
            kernel = fspecial('log');
    end
    
    if filtertype > 0
        % Build shader from kernel:
        shader = EXPCreateStatic2DConvolutionShader(kernel, 3, 3, 0, 1);
    else
        % No filtering requested: Select fixed-function pipeline
        shader = 0;
    end
    
    % Open video capture device:
    grabber = Screen('OpenVideoCapture', win, deviceIndex);
    
    % Start video capture:
    Screen('StartVideoCapture', grabber);
    
    count = 0;
    t=GetSecs;
    
    % Repeat until keypress or timeout of 10 minutes:
    while ((GetSecs - t) < 600) && ~KbCheck
        
        % Wait for next video image:
        tex = Screen('GetCapturedImage', win, grabber, 1);
        
        % Valid?
        if tex > 0
            % Draw new video texture from framegrabber.
            Screen('DrawTexture', win, tex, [], [], [], [], [], [], shader);
            Screen('Close', tex);
            
            % Show it:
            Screen('Flip', win);
        end
        
        count = count + 1;
    end
    
    telapsed = GetSecs - t
    
    % Stop and close video source:
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
    
    % Close down everything else:
    Screen('CloseAll');
    
    avgfps = count / telapsed
catch
    sca;
    psychrethrow(psychlasterror);
end
