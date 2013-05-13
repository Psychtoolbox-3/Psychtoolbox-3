function VideoOfflineCaptureDemo(fullscreen, fullsize, roi, deviceId)
% Demonstrate use of built-in GStreamer video capture engine to capture first
% into memory, then retrieve corresponding video textures after end of capture.
%
% VideoOfflineCaptureDemo([fullscreen=0][, fullsize=0][, roi=[0 0 640 480]][,deviceId=0])
%
% VideoOfflineCaptureDemo initializes the first attached and supported camera on
% your computer (e.g, the built-in iSight of Apple Macintosh computers),
% then records video from it into memory until you press a key. It then
% fetches and displays all previously recorded images, then quits. Abort
% display by pressing any key on the keyboard.
%
% By default, the maximum supported framerate is requested and the
% timecode and interframe interval of each captured image is displayed in
% the top-left corner of the display.
%
% See also ImagingVideoCaptureDemo, VideoDelayloopMiniDemo and a few other
% nice demos.
%
% Optional parameters:
%
% 'fullscreen' If set to non-zero value, the image is displayed in a
% fullscreen window, as usual, otherwise a normal GUI window is used.
%
% 'fullsize' If set to 1, the cameras image is scaled up to full screen
% resolution, ie. so it fills the maximum amount of display area, but
% preserving the original aspect ratio.
%
% 'roi' Selects a rectangular subregion of the camera for display. By
% default, it selects a [0 0 640 480] rectangle, ie. the full are of a
% camera with 640 x 480 pixels resolution. This parameter may need tweaking
% for some cameras, as some drivers have bugs and don't work well with all
% settings.
%
% 'deviceId' Device index of video capture device. Defaults to system default.
%

% History:
% 17.04.2011  mk   Written.

AssertOpenGL;

if nargin < 1
    fullscreen=[];
end

if isempty(fullscreen)
    fullscreen=0;
end;

if nargin < 2
    fullsize=[];
end

if isempty(fullsize)
    fullsize=0;
end

if nargin < 3
    roi = [];
end

if nargin < 4 || isempty(deviceId)
    deviceId = [];
end

screenid=max(Screen('Screens'));

try
    if fullscreen<1
        win=Screen('OpenWindow', screenid, 0, [0 0 800 600]);
    else
        win=Screen('OpenWindow', screenid, 0);
    end;

    % Initial flip to a blank screen:
    Screen('Flip',win);

    % Set text size for info text. 24 pixels is also good for Linux.
    Screen('TextSize', win, 24);
    
    % Open framegrabber on specified capture device 'deviceId' with specified
    % resolution and region of interest 'roi'. We enforce use of capture engine
    % type '3', the GStreamer engine, as only this engine supports offline capture:
    grabber = Screen('OpenVideoCapture', win, deviceId, roi, [], [], [], [], [], 3);

    % Start videocapture at maximum framerate, do not drop frames, but keep'em all:
    [fps, t] = Screen('StartVideoCapture', grabber, realmax, 0);

    % Wait for user to stop capture:
    DrawFormattedText(win, 'Video capture active.\nPress any key to stop capture\nand start display of the results.', 'center', 'center', 255);
    Screen('Flip', win);
    KbStrokeWait;

    % Stop video capture without discarding all internally stored frames, so we can fetch'em later:
    Screen('StopVideoCapture', grabber, 0);

    % Fetch loop: Fetches and displays all captured images until keypress or all images
    % are used up:
    dstRect = [];
    oldpts = 0;
    count = 0;
    tex=0;

    while 1
        if KbCheck
            break;
        end;
        
        [tex pts nrpending] = Screen('GetCapturedImage', win, grabber, [], tex); %#ok<NASGU>
        
        if tex > 0
	    % fprintf('tex = %i  pts = %f nrpending = %i\n', tex, pts, nrpending);

            % Perform first-time setup of transformations, if needed:
            if fullsize && (count == 0)
                texrect = Screen('Rect', tex);
                winrect = Screen('Rect', win);
                sf = min([RectWidth(winrect) / RectWidth(texrect) , RectHeight(winrect) / RectHeight(texrect)]);
                dstRect = CenterRect(ScaleRect(texrect, sf, sf) , winrect);
            end

            % Draw new texture from framegrabber.
            Screen('DrawTexture', win, tex, [], dstRect);

            % Print pts as delta since start of capture:
            Screen('DrawText', win, sprintf('%.4f', pts - t), 0, 0, 255);
            if count>0
                % Compute delta:
                delta = (pts - oldpts) * 1000;
                oldpts = pts;
                % Display delta between successive frames in msecs:
                Screen('DrawText', win, sprintf('%.4f', delta), 0, 20, 255);
            end;
            
            % Show it.
            Screen('Flip', win);
	else
            % Done with fetching.
	    break;
        end;

        count = count + 1;
    end;

    % Close capture device: This will release all not yet fetched images:
    Screen('CloseVideoCapture', grabber);

    % Close window, clean up.
    Screen('CloseAll');
catch
    Screen('CloseAll');
end;
