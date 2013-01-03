function KinectDemo
% KinectDemo - Capture and display video and depths data from a Kinect box.
%
% This connects to a Microsoft Kinect device on the USB bus, then captures
% and displays video and depths data delivered by the Kinect.
%
% Press any key to quit the demo.
%
% See help PsychKinect, help PsychKinectCore and help InstallKinect for
% further info.
%

% History:
% 25.11.2010  mk  Written.

AssertOpenGL
%PsychDebugWindowConfiguration;
InitializeMatlabOpenGL([],[],1);
GetSecs;

dst1 = [0, 0, 640, 480];
dst2 = [650, 0, 650+640, 480];

w = Screen('Openwindow', 0, 0, [0 0 1300 500], [], [], [], [], [], kPsychGUIWindow);
kinect = PsychKinect('Open');
PsychKinect('Start', kinect);
count = 0;
ts = GetSecs;

while 1
    [rc, cts] = PsychKinect('GrabFrame', kinect, 1);
    if rc > 0
        count = count + 1;
        %fprintf('Kinect frame %i, cts = %f\n', count, cts);
        [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 0, 1);
        if width > 0 && height > 0
            tex = Screen('SetOpenGLTextureFromMemPointer', w, [], imbuff, width, height, channels, 1, GL.TEXTURE_RECTANGLE_EXT);
            Screen('DrawTexture', w, tex, [], dst1);
            Screen('Close', tex);
        end
        
        [imbuff, width, height, channels] = PsychKinect('GetImage', kinect, 1, 1);
        if width > 0 && height > 0
            tex = Screen('SetOpenGLTextureFromMemPointer', w, [], imbuff, width, height, channels, 1, GL.TEXTURE_RECTANGLE_EXT);
            Screen('DrawTexture', w, tex, [], dst2);
            Screen('Close', tex);
        end
        
        if 0
            [imbuff, width, height, channels] = PsychKinect('GetDepthImage', kinect, 2, 0);
            if 0
                imagesc(imbuff);
                drawnow;
                minz = min(min(imbuff))
                maxz = max(max(imbuff))
            else
                scatter3(imbuff(1,:,:), imbuff(2,:,:), imbuff(3,:,:));
                drawnow;
            end
        end
        
        PsychKinect('ReleaseFrame', kinect);
        Screen('Flip', w);
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

PsychKinect('Stop', kinect);
PsychKinect('Close', kinect);
sca;
