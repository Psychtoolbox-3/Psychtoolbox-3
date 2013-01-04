function KinectRawRecordingDemo(deviceId)
% KinectRawRecordingDemo - Capture, display and record video and depths data from a Kinect box.
%
% This connects to a Microsoft Kinect device on the USB bus, then captures
% and displays video and depths data delivered by the Kinect. It also stores
% the data to filesystem
% Press any key to quit the demo.
%
% See help PsychKinect, help PsychKinectCore and help InstallKinect for
% further info.
%

% History:
% 25.11.2010  mk  Written.

if nargin < 1
    deviceId = [];
end

AssertOpenGL
%PsychDebugWindowConfiguration;
InitializeMatlabOpenGL([],[],1);
GetSecs;

dst1 = [0, 0, 640, 480];
dst2 = [650, 0, 650+640, 480];

w = Screen('Openwindow', 0, 0, [0 0 1300 500]);

% Open 'deviceId'th connected Kinect.
% Use 300 internal buffers for buffering up to 10 secs of
% capture at 30 fps max in the fifo. Disable (0) bayer filtering
% to save memory:
kinect = PsychKinect('Open', deviceId, 300, 0);

PsychKinect('Start', kinect);
count = 0;
ts = GetSecs;
elapsed = [];

while 1
    [rc, cts] = PsychKinect('GrabFrame', kinect, 1);
    if rc > 0
        count = count + 1;

        % Visualization:

        % Image:
        if 0
            %fprintf('Kinect frame %i, cts = %f\n', count, cts);
            [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 0, 1);
            if width > 0 && height > 0
                tex = Screen('SetOpenGLTextureFromMemPointer', w, [], imbuff, width, height, channels, 1, GL.TEXTURE_RECTANGLE_EXT);
                Screen('DrawTexture', w, tex, [], dst1);
                Screen('Close', tex);
            end
        end

        % Color coded depth map:
        if 0
            [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 1, 1);
            if width > 0 && height > 0
                tex = Screen('SetOpenGLTextureFromMemPointer', w, [], imbuff, width, height, channels, 1, GL.TEXTURE_RECTANGLE_EXT);
                Screen('DrawTexture', w, tex, [], dst2);
                Screen('Close', tex);
            end
        end

        Screen('Flip', w, [], 2, 2);

        % Image & Depth data disc writeout:
        if 1
            t1 = GetSecs;
            [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 0);
            depth = PsychKinect('GetDepthImage', kinect, 8);
            if ~IsWin
                save(sprintf('/tmp/testdepths_%i.mat', count), 'depth', 'imbuff');
            else
                save(sprintf('D:\\temp\\testdepths_%i.mat', count), 'depth', 'imbuff');
            end
            elapsed(end+1) = (GetSecs - t1) * 1000;
        end

        PsychKinect('ReleaseFrame', kinect);
    else
        WaitSecs('YieldSecs', 0.005);
    end

    [x,y,buttons]=GetMouse;
    if (x == 0 && y == 0)
        GetClicks;
    end

    if KbCheck(-1)
        break;
    end
end

fprintf('Average fps = %f [%i]\n', count / (GetSecs - ts), count);

twrite = mean(elapsed)
close all ; plot(elapsed);
PsychKinect('Stop', kinect);
PsychKinect('Close', kinect);
sca;

% load /tmp/testdepths_xxx.mat
% imagesc(transpose(double(depth)))
% img = transpose(imbuff); imshow(img)
