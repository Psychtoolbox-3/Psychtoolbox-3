function VideoCaptureLatencyTest(texfetch, fullscreen, depth)
% VideoCaptureLatencyTest
%
% This routine tests the latency of Videocapture for
% a given Camera -> Computer -> PTB -> gfx-card -> display combo.
%
% Working principle and usage:
%
% 1. Point camera to monitor so that cameras field of view captures monitor
% display area.
%
% 2. Calibration: First a black screen is shown to the cam and captured to
% find the average intensity of an image of a black display.
% Then (over ten trials), the screen turns white and the script measures
% the time delay between onset of white screen and arrival of a captured
% image with at least 25% more intensity than the 'black'-image -- our
% white image. This procedure is repeated 10 times and latencies averaged.
%
% -> Rough estimate of delay from scene change to arrival of frame.
% -> You have to add the display onset delay if captured image should be
% shown to subject for video feedback loops.
%
% Parameters:
%
% texfetch == 0 Measure only capture latency. 1 == Measure texture fetch
% latency as well. 2 == Measure drawing and onset latency as well.
%
% History:
% 2/9/06 mk Written.

% Running OpenGL-PTB? This is required.
AssertOpenGL;

rate = 30;

% Find stimulus display:
screen=max(Screen('Screens'));

if nargin < 1
    texfetch=0;
end;

texfetch

if nargin < 2
   fullscreen=1;
end;

fullscreen

if nargin < 3
   depth = 1;
end;
depth

if IsOctave
    more off;
    ignore_function_time_stamp('all');
end

%try
    % Open display with default settings:
    if fullscreen<1
	    win=Screen('OpenWindow', screen, 0, [0 0 800 600]);
    else
       win=Screen('OpenWindow', screen, 0);
    end;
    
    [w h]=Screen('WindowSize', win);
    
    % Setup font and textsize: 24 is also a good size for Linux.
    Screen('TextSize', win, 24);
    Screen('TextStyle', win, 1);
    
    % Show initial blank screen:
    Screen('Flip',win);
    
    % Open default video capture device:
    grabber = Screen('OpenVideoCapture', win, 0, [0 0 640 480], depth);
    
    % Start video capture: We request at least a 30 Hz capture rate and lowlatency (=1) mode:
    % Requested start time is 4 seconds from now:
    [fps capturestarttime]=Screen('StartVideoCapture', grabber, rate, 1, GetSecs + 4);
    fps
    
    oldpts = 0;
    count = 0;
    t=GetSecs;
    while (GetSecs - t) < 600
        if KbCheck
            break;
        end;
        
        [tex pts nrdropped intensity]=Screen('GetCapturedImage', win, grabber, 0);
        % fprintf('tex = %i  pts = %f\n', tex, pts);
        if (tex>0)
            % Print instructions:
            Screen('DrawText', win, 'Point camera to center of Display - Then press any key', 10, h-50);
            % Print pts:
            Screen('DrawText', win, ['Capture timestamp (s): ' sprintf('%.4f', pts)], 0, 0, 255);
            if count>0
                % Compute delta:
                delta = (pts - oldpts) * 1000;
                oldpts = pts;
                Screen('DrawText', win, ['Delta since last frame (ms): ' sprintf('%.4f', delta)], 0, 24, 255);
                Screen('DrawText', win, ['Delta start to first frame (ms): ' sprintf('%.4f', camstartlatency)], 0, 96, 255);
            else
                % First image in video fifo.
                camstartlatency = (pts - capturestarttime) * 1000.0
            end;
            
            % Compute and print intensity:
            Screen('DrawText', win, ['Mean intensity: ' sprintf('%.2f', intensity)], 0, 48, 255);
           
            % Print count of dropped frames since last fetch:
            Screen('DrawText', win, ['Dropped frames: ' sprintf('%i', nrdropped)], 0, 72, 255);

            % New texture from framegrabber.
            Screen('DrawTexture', win, tex);
            Screen('Flip', win, 0, 0, 2);
            Screen('Close', tex);
            tex=0;
            count = count + 1;
        end;        
    end;
    telapsed = GetSecs - t
    
    calcedlatency = PsychCamSettings('EstimateLatency', grabber) * 1000.0
    tminimum = 10000000000.0;
    tminimum2 = 10000000000.0;
    tcount = 1;


    for reps=1:5
    % Ok, preview finished. Cam should now point to the screen and
    % we can start the latency measurement procedure:
    ntrials = 10;
    tavgdelay = 0;
    threshold = 0;
    
    for i=1:ntrials
        % Make screen black:
        Screen('FillRect', win, 0);
        Screen('Flip', win);
        
        % Capture at least 20 frames to make sure we're really black...
        for j=1:20
           Screen('GetCapturedImage', win, grabber, 2);
        end;
       
        % Capture a black frame:
        tex=0;
        while (tex==0)
            [tex pts nrdropped intensity]=Screen('GetCapturedImage', win, grabber, 0);
        end;

        if (tex>0) Screen('Close', tex); end;
    
%        if threshold==0
          blackintensity = intensity;
          threshold = blackintensity * 1.02;
%        end;
     
        % Eat up all enqueued images, if any.
        tex=1;
        while (tex>0)
            [tex newpts]=Screen('GetCapturedImage', win, grabber, 0);
            if (tex>0)
                Screen('Close', tex);
                pts = newpts;
            end;
        end;

        % Make display white:
        Screen('FillRect', win, 255);
        % Show it and take onset timestamp:
        tonset = Screen('Flip', win);
        
        % Now we loop and capture until intensity exceeds threshold:
        tex=0;
        intensity=0;
        lagframes=0;
        ptsblack = pts;
        
        while (tex>=0 && intensity < threshold)
            if texfetch>0
                waitmode=1;
            else
                waitmode=2;
            end;
            
            [tex pts nrdropped intensity]=Screen('GetCapturedImage', win, grabber, waitmode);
            
            % Captured! Take timestamp:
            tdelay=(GetSecs - tonset) * 1000;
        
            if tex>0
               lagframes=lagframes+1;
               if texfetch>1 && intensity>=threshold
                  Screen('DrawTexture',win,tex,[],[],[],0);
                  tdelay = (Screen('Flip',win) - tonset) * 1000;
               end;
               Screen('Close', tex);
               tex=0;
            end;
        end;
                
        tdelay
        lagframes
        ptslag = (pts - ptsblack) * 1000
        
        tavgdelay=tavgdelay + tdelay;
        tminimum2=min(tminimum2, tdelay);      
        ttotal(tcount)=tdelay;
        tcount = tcount + 1;

        % Next trial...
    end;
        
    % Stop and restart capture loop:
    Screen('StopVideoCapture', grabber);

    % Restart with a random delay between 0.0 and 1.0 seconds...
    Screen('StartVideoCapture', grabber, rate, 1, GetSecs + rand);

    tavgdelay = tavgdelay / ntrials;
    tavgdelay

    tminimum = min(tminimum, tavgdelay);
    tdelays(reps) = tavgdelay;

    end

    % Shutdown capture loop, close screens...
    Screen('CloseVideoCapture', grabber);
    Screen('CloseAll');

    tminimum
    tminimum2

    figure;
    plot(tdelays);
    figure;
    plot(ttotal);

    return;
    
    %catch
   %Screen('CloseAll');
   %psychrethrow(psychlasterror);
%end;
