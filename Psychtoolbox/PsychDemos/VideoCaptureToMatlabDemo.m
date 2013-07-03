function VideoCaptureToMatlabDemo(deviceIndex)
% VideoCaptureToMatlabDemo([deviceIndex])
%
% Minimalistic demo on how to capture video data and return it in a
% Matlab/Octave matrix:

AssertOpenGL;

if nargin < 1
    deviceIndex = [];
end

try
    % Open a minimalistic window (only 10 x 10 pixels). This is needed
    % because we need to provide a 'win'dowhandle to the videocapture
    % functions:
    win=Screen('OpenWindow', 0, 0, [0 0 10 10]);

    % Open videocapture device:
    grabber = Screen('OpenVideoCapture', win, deviceIndex, [0 0 640 480]);
    % Start capture with requested 30 fps:
    Screen('StartVideoCapture', grabber, 30, 1);

    oldpts = 0;
    count = 0;
    t=GetSecs;
    % Capture for 600 seconds or until keypress:
    while (GetSecs - t) < 600 
        if KbCheck
            break;
        end;
        
        % Retrieve next captured image in 'rawImage'. The 'waitforImage=2'
        % flag disables texture creation, so 'tex' is actually an empty
        % handle. The 'specialMode'=2 flag requests video data as matrix:
        [tex pts nrdropped rawImage]=Screen('GetCapturedImage', win, grabber, 2, [], 2);

        % Bits of accounting and stats to the Matlab window:
        if count>0
            % Compute delta:
            delta = (pts - oldpts) * 1000;
            oldpts = pts;
            fprintf('%.4f\n', delta);
        end;

        % Convert rawImage matrix into a matrix suitable for display with
        % Matlabs imshow(). imshow needs a height x width x 3 colors
        % matrix, whereas rawImage is a c by width x height matrix with c=1
        % for luminance data, c=3 or 4 for RGB data, where the 4th
        % component - if present - is a useless alpha channel.
        channels = min(size(rawImage,1), 3);
        for ci=1:channels
            if channels == 1
                tci=1;
            else
                tci = 4 - ci;
            end
            matImage(1:size(rawImage,3), 1:size(rawImage,2), tci) = transpose(squeeze(rawImage(ci,:,:)));
        end
        
        % Show image in a figure window:
        imshow(matImage);
        drawnow;
        count = count + 1;
    end;
    
    telapsed = GetSecs - t
    
    % Stop capture, close engine and onscreen window:
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
    Screen('CloseAll');
    avgfps = count / telapsed
catch
   Screen('CloseAll');
end;
