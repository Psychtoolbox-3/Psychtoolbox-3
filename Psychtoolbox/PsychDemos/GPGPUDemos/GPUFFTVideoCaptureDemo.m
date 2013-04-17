function GPUFFTVideoCaptureDemo(usegpu, showfft, fwidth, roi, depth, deviceId, cameraname)
% GPUFFTVideoCaptureDemo - Demonstrate use of GPGPU computing for live filtering via 2D-FFT.
%
% This demo makes use of the FOSS GPUmat toolbox to perform a GPU
% accelerated 2D FFT + filtering in frequency space + 2D inverse FFT on a
% live video feed from video capture or movie playback. GPUmat allows to
% use NVidia's CUDA gpu computing framework on supported NVidia gpu's
% (GeForce-8000 series and later, aka Direct3D-10 or OpenGL-3 capable).
%
% It shows how a Psychtoolbox floating point texture (with video content
% inside) can be efficiently passed to GPUmat as a matrix of GPUsingle data
% type, which is stored and processed on the GPU. Then it uses GPUmat's fft
% routines for forward/inverse fft's and matrix manipulation. Then it
% returns the final image to Psychtoolbox as a new floating point texture
% for display. The same steps are carried out with Matlab/Octave's regular
% fft routines on the cpu for reference.
%
% Requires the freely downloadable NVidia CUDA-5.0 SDK/Runtime and the free
% and open-source GPUmat toolbox as well as a compute capable NVidia
% graphics card.
%
% Usage:
%
% GPUFFTVideoCaptureDemo([usegpu=1][, showfft=0][, fwidth=11][, roi=[0 0 640 480]][, depth=1][, deviceId=0][, cameraname])
%
% Parameters:
%
% 'usegpu' = 0 to use regular Matlab/Octave fft on CPU, 1 = to use GPUmat on GPU.
%
% 'showfft' = 1 to show amplitude spectrum of video in usegpu=1 mode.
%
% 'fwidth' = Width of low-pass filter kernel in frequency space units.
%
% 'roi' Selects a rectangular subregion of the camera for display. By
% default, it selects a [0 0 640 480] rectangle, ie. the full area of a
% camera with 640 x 480 pixels resolution. This parameter may need tweaking
% for some cameras, as some drivers have bugs and don't work well with all
% settings.
%
% 'depth' = 1 for Mono, 3 for color processing.
%
% 'deviceId' Device index of video capture device. Defaults to system default.
%
% 'cameraname' Name string for selection of video capture device. This is
% only honored if 'deviceId' is a negative number, and only for certain
% video capture plugins. Defaults to none.
%

% History:
% 5.02.2013  mk  Written.

showmask = 0;

if nargin < 1 || isempty(usegpu)
    usegpu = 1;
end

if nargin < 2 || isempty(showfft)
    showfft = 0;
end

% Default filter width is 11 units:
if nargin < 3 || isempty(fwidth)
    fwidth = 11;
end

if nargin < 4
    roi = [];
end

if nargin < 5 || isempty(depth)
    depth = 1;
end

if nargin < 6
    deviceId = [];
end

if nargin < 7
    cameraname = [];
end

% At least on Linux + NVidia + CUDA-5.0 it is currently safe to keep
% texture resources mapped during OpenGL operations on CUDA-OpenGL interop
% textures. On Windows we don't know yet. On OSX it is a total no-go. However,
% the CUDA-5.0 programming guide strongly advices against keeping resources
% mapped during OpenGL operations, as the "results of such actions will be
% undefined". In other words, nothing guarantees it will work, or even if it
% works now, will continue to work in the future. So better safe than sorry
% and default this to off on all platforms. Adventureous or performance
% desparate users can enable it on their system and see if it works.
keepmapped = 0;

% Running under PTB-3 hopefully?
AssertOpenGL;

% Skip timing tests and calibrations for this demo:
oldskip = Screen('Preference','SkipSyncTests', 2);

% Open onscreen window with black background on (external) screen,
% enable GPGPU computing support:
screenid = max(Screen('Screens'));

PsychImaging('PrepareConfiguration');
% We explicitely request the GPUmat based api:
PsychImaging('AddTask', 'General', 'UseGPGPUCompute', 'GPUmat');
win = PsychImaging('OpenWindow', screenid, 0);

try
    Screen('TextSize', win, 28);
    
    % Color or mono capture and processing?
    if ~usegpu
        depth = 1;
    end
    
    if depth > 2
        mc = 3;
    else
        mc = 1;
    end
    
    % Open videocapture device:
    grabber = Screen('OpenVideoCapture', win, deviceId, roi, depth, [], [], cameraname);
    Screen('StartVideoCapture', grabber, realmax, 1);
    
    tex = Screen('GetCapturedImage', win, grabber, 1);
    
    % Perform first-time setup of transformations on first frame:
    % Get size of a video frame:
    texrect = Screen('Rect', tex);
    
    % Get window size:
    winrect = Screen('Rect', win); 
    
    if showfft
        dstRect = texrect;
        dstRect = OffsetRect(dstRect, 40, 40);
        fftRect = AdjoinRect(texrect, dstRect, RectRight);
        fftRect = OffsetRect(fftRect, 40, 40);
    else
        % Compute scaling factor and dstRect to upscale video to fullscreen display:
        sf = min([RectWidth(winrect) / RectWidth(texrect) , RectHeight(winrect) / RectHeight(texrect)]);
        dstRect = CenterRect(ScaleRect(texrect, sf, sf) , winrect);
    end
    
    % Define initial "amplitude spectrum multiplicative mask"
    % for storage of the frequency filter function in the
    % frequency domain (post FFT, pre IFFT):
    hh = RectHeight(texrect) / 2;
    hw = RectWidth(texrect) / 2;
    
    % fwidth controls frequency - lower values = lower cutoff frequency:
    [x,y] = meshgrid(-hw:hw,-hh:hh);
    mask = exp(-((x/fwidth).^2)-((y/fwidth).^2));
    
    % Quick and dirty trick: Cut off a bit from mask, so size of mask and
    % frequency spectrum matches -- Don't do this at home (or for real research scripts)!
    mask = mask(2:end, 2:end);
    
    % Show it in figure:
    if showmask
        figure; %#ok<*UNRCH>
        imshow(mask);
        title('Gaussian low pass filter mask for FFT space filtering:');
    end
    
    % Apply inverse fftshift on mask, so we don't need to
    % fftshift/ifftshift the FFT video input data. Instead we
    % shift the mask into the format/position of the
    % input/output post-FFT/pre-IFFT stream.
    mask = ifftshift(mask);
    
    if showmask
        figure;
        imshow(mask);
        title('ifftshifted Gaussian low pass filter mask for FFT space filtering:');
    end
    
    if usegpu
        % Turn our filter 'mask' from cpu version above into a matrix on GPU:
        % We must also transpose the mask, because 'mask' was generated to
        % match the shape of the fs matrix on the cpu in Matlab/Octave. As
        % Matlab and Octave use column-major storage, whereas Psychtoolbox uses
        % row-major storage, all the GPU variables we got from Psychtoolbox
        % textures are transposed with respect to the Matlab version. A simple
        % transpose fixes this for our mask to match GPU format:
        FM = transpose(GPUsingle(mask));
        
        % Create 32 bpc floating point resolution offscreen
        % window as intermediate storage:
        texf = Screen('OpenOffscreenWindow', win, 0, texrect, 128);

        if mc > 1
          R = ones(4, size(FM, 1), size(FM, 2), GPUsingle);
        else
          R = ones(1, size(FM, 1), size(FM, 2), GPUsingle);
        end
    end
    
    % First texture was only for setup above. Get rid of it:
    Screen('Close', tex);

    tex = [];
    A = [];
    rtex = [];
    fftmag = [];
    count = 0;
    
    % Take start timestamp for benchmarking:
    tstart = GetSecs;
    
    % Video capture + processing + display loop: Runs until key press.
    while ~KbCheck        
        %% Perform the FFT -> multiply -> IFFT filtering in frequency domain:
        if usegpu
            %% Use GPUmat/CUDA for GPU accelerated processing:

            % Wait for next video image to arrive and retrieve a 'tex'ture
            % handle to it. Recylcle 'tex'ture from previous cycle, if any:
            tex = Screen('GetCapturedImage', win, grabber, 1, tex);
            
            % Our video 'tex' is in 8 bit integer format, but GPUmat
            % needs single() format aka 32 bpc float. We convert the
            % image by drawing it one-to-one into texf, a 32 bpc float
            % RGBA offscreen window:
            Screen('DrawTexture', texf, tex, [], [], [], 0);
            
            % Create GPUsingle matrix with input image from float tex with video image:
            A = GPUTypeFromToGL(0, texf, [], A, keepmapped);
            
            for c = 1:mc
                % Extract 1st layer, the red aka luminance channel from it for further use:
                % Note: The 1st dimension of a GPUsingle matrix created from a
                % Psychtoolbox texture always selects the "color channel":
                AL = A(c, :, :);
                
                % Squeeze it to remove the singleton 1st dimension of A, because
                % fft2 and ifft2 can only work on 2D input matrices, not 3D
                % matrices, not even ones with a singleton dimension, ie., a 2D
                % matrix in disguise:
                AL = squeeze(AL);
                
                % Perform forward 2D FFT on GPU:
                F = fft2(AL);
                
                if showfft
                    % Generate the amplitude spectrum, scaled to 1/100th via abs(FS)/100.0,
                    % then converte the image into a texture 'fftmag' and display it:
                    fftmag = GPUTypeFromToGL(1, abs(fftshift(F))/100.0, [], fftmag, keepmapped);
                    Screen('Blendfunction', win, [], [], [double(c==1), double(c==2), double(c==3), 1]);
                    Screen('DrawTexture', win, fftmag, [], fftRect, [], 0);
                    Screen('Blendfunction', win, [], [], [1 1 1 1]);
                    DrawFormattedText(win, 'Amplitude spectrum of video:', fftRect(RectLeft), fftRect(RectTop) - 30, [0 255 0]);
                end
                
                % Filter the amplitude spectrum by point-wise multiply with filter FM:
                F = F .* FM;
                
                % Process inverse 2D FFT on GPU:
                B = ifft2(F);
                
                % Extract real component for display:
                R(c,:,:) = real(B);
            end
            
            % Convert R back into a floating point luminance texture 'tex' for
            % processing/display by Screen. Here it is fine to pass a
            % single-layer matrix for getting a luminance texture, because we
            % don't use Screen('MakeTexture') or similar and don't perform the
            % normalization step which would convert into a troublesome RGB
            % texture.
            rtex = GPUTypeFromToGL(1, R, [], rtex, keepmapped);
            
            % Show final image onscreen, scaled to full window size:
            Screen('DrawTexture', win, rtex, [], dstRect);
            if showfft
                DrawFormattedText(win, 'GPU filtering', dstRect(RectLeft), dstRect(RectTop) - 30, [0 255 0]);
            else
                DrawFormattedText(win, 'GPU filtering', 'center', 30, [0 255 0]);
            end
        else
            %% Use classic Matlab/Octave fft path on cpu:
            
            % Wait for next video image to arrive and retrieve a 'tex'ture
            % handle to it. Recylcle 'tex'ture from previous cycle, if any:
            %tex = Screen('GetCapturedImage', win, grabber, 1, tex);

            % Retrieve next captured image in 'rawImage'. The 'waitforImage=2'
            % flag disables texture creation, so 'tex' is actually an empty
            % handle. The 'specialMode'=2 flag requests video data as matrix:
            [tex pts nrdropped rawImage] = Screen('GetCapturedImage', win, grabber, 2, [], 2); %#ok<ASGLU>
            
            % Convert uint8 luminance image 'rawImage' into Matlab 2D
            % single() matrix with luminance data:
            mframe = single(transpose(squeeze(rawImage))) / 255;
            
            % Perform 2D-FFT in Matlab/Octave on CPU:
            f = fft2(mframe);
            
            % Low pass filter by point-wise multiply of fs with filter 'mask' in
            % frequency space:
            f = f .* mask;
            
            % Perform inverse 2D-FFT in Matlab/Octave on CPU:
            b = ifft2(f);
            
            % Extract real component for display:
            r = real(b);

            % Convert real component result image into texture:
            
            % Doing it as float texture is a little bit slower:
            %rtex = Screen('MakeTexture', win, double(r), [], [], 2);
            
            % Doing it as uint8 texture is a tiny bit faster:
            rtex = Screen('MakeTexture', win, uint8(r * 255));
            
            % Show final image onscreen, scaled to full window size:
            Screen('DrawTexture', win, rtex, [], dstRect);
            DrawFormattedText(win, 'CPU filtering', 'center', 30, [255 255 0]);
            Screen('Close', rtex);
        end
        
        % Show it at next vertical retrace, don't clear framebuffer or
        % do anything else:
        Screen('Flip', win, [], 2, 2);
        
        % Increment frame counter:
        count = count + 1;
    end % Next capture loop iteration.
    
    % Take end timestamp for benchmarking:
    tend = GetSecs;
    fps = count / (tend - tstart);
    fprintf('Average FPS %f\n', fps);
    
    % Reset interop cache before closing the windows and textures:
    GPUTypeFromToGL(4);
    
    % Stop and close video capture engine:
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
    
    % Close window, which will also release all textures:
    sca;
    
    % Restore sync test setting:
    Screen('Preference','SkipSyncTests', oldskip);
    
catch %#ok<CTCH>
    % Error handling.
    
    % Reset interop cache before closing the windows and textures:
    GPUTypeFromToGL(4);
    
    % Close window, this also shuts down video capture:
    sca;
    
    % Restore sync test setting:
    Screen('Preference','SkipSyncTests', oldskip);
    
    % Propagate the error upwards:
    psychrethrow(psychlasterror);
end

% Done.
end
