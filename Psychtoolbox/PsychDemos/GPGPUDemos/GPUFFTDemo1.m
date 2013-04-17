function GPUFFTDemo1(fwidth)
% GPUFFTDemo1([fwidth=11]) - Demonstrate use of GPGPU computing for 2D-FFT.
%
% This demo makes use of the FOSS GPUmat toolbox to perform a GPU
% accelerated 2D FFT + filtering in frequency space + 2D inverse FFT on our
% favourite bunnies. GPUmat allows to use NVidia's CUDA gpu computing
% framework on supported NVidia gpu's (GeForce-8000 series and later, aka
% Direct3D-10 or OpenGL-3 capable).
%
% It shows how a Psychtoolbox floating point texture (with our bunnies
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
% Parameters:
%
% 'fwidth' = Width of low-pass filter kernel in frequency space units.
%

% History:
% 5.02.2013  mk  Written.

% Default filter width is 11 units:
if nargin < 1 || isempty(fwidth)
    fwidth = 11;
end

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
w = PsychImaging('OpenWindow', screenid, 0);

try
    Screen('TextSize', w, 28);
    DrawFormattedText(w, 'Please be patient throughout the demo.\nSome benchmark steps are time intense...', 'center', 'center', 255);
    Screen('Flip', w);
    
    % Read our beloved bunny image from filesystem:
    bunnyimg = imread([PsychtoolboxRoot 'PsychDemos/konijntjes1024x768.jpg']);
    
    % Add a fourth neutral alpha channel, so we can create a 4-channel
    % RGBA texture. Why? Some GPU compute api's, e.g., NVidia's CUDA,
    % don't like RGB textures, so we need to do this to avoid errors:
    bunnyimg(:,:,4) = 255;
    
    % Turn it into a double matrix for conversion into a floating point
    % texture, and normalize/convert its color range from 0-255 to 0-1:
    dbunny = double(bunnyimg)/255;
    
    % Maketexture a texture out of it in float precision with upright
    % texture orientation, so further processing in teamwork between
    % GPUmat and Psychtoolbox is simplified. Note: This conversion
    % would turn a luminance texture into a troublesome RGB texture,
    % which would only work on Linux and maybe Windows, but not on OSX.
    % We avoid this by never passing in a luminance texture if it is
    % intended to be used with GPUmat aka NVidia's CUDA framework.
    bunnytex = Screen('MakeTexture', w, dbunny, [], [], 2, 1);
    
    %% Do the FFT -> manipulate -> IFFT cycle once in Matlab/Octave as a
    % reference:
    
    % Extract 2nd layer, the green color channel as an approximation of
    % luminance:
    dbunny = dbunny(:,:,2);
    
    % Show it pre-transform:
    close all;
    imshow(dbunny);
    title('Pre-FFT Bunny.');
    
    % Perform 2D-FFT in Matlab/Octave on CPU:
    f = fft2(dbunny);

    % Shift DC component to center of spectrum "image":
    fs = fftshift(f);

    % Define a low-pass filter in frequency space, ie., an amplitude mask
    % to point-wise multiply with the FFT spectrum of FFT transformed
    % image:
    % fwidth controls frequency - lower values = lower cutoff frequency:
    hh = size(fs, 1)/2;
    hw = size(fs, 2)/2;
    [x,y] = meshgrid(-hw:hw,-hh:hh);
    mask = exp(-((x/fwidth).^2)-((y/fwidth).^2));
    
    % Quick and dirty trick: Cut off a bit from mask, so size of mask and
    % frequency spectrum matches -- Don't do this at home (or for real research scripts)!
    mask = mask(2:end, 2:end);
    
    figure;
    imshow(mask);
    title('Gaussian low pass filter mask for FFT space filtering:');

    tcpu = GetSecs;
    
    for trials = 1:10
        % Low pass filter by point-wise multiply of fs with filter 'mask' in
        % frequency space:
        fs = fs .* mask;
        
        % Invert shift of filtered fs:
        f = ifftshift(fs);
        
        % Perform inverse 2D-FFT in Matlab/Octave on CPU:
        b = ifft2(f);
    end
    
    tcpu = (GetSecs - tcpu) / trials;
    fprintf('Measured per trial runtime of CPU FFT: %f msecs [n=%i trials].\n', tcpu * 1000, trials);
    
    % Extract real component for display:
    r = real(b);
    
    figure;
    imshow(r);
    title('Post-FFT->Filter->IFFT Bunny.');
    
    %% Perform FFT->Process->IFFT on GPU via GPUmat / CUDA:
    
    % First Show input image:
    Screen('DrawTexture', w, bunnytex);
    Screen('TextSize', w, 28);
    DrawFormattedText(w, 'Original pre-FFT Bunny:\nPress key to continue.', 'center', 40, [255 255 0]);
    Screen('Flip', w);
    KbStrokeWait(-1);
    
    % Create GPUsingle matrix with input image from RGBA float bunnytex:
    A = GPUTypeFromToGL(0, bunnytex, [], [], 0);
    
    % Extract 2nd layer, the green channel from it for further use:
    % Note: The 1st dimension of a GPUsingle matrix created from a
    % Psychtoolbox texture always contains the "color channel":
    A = A(2, :, :);
    
    % Squeeze it to remove the singleton 1st dimension of A, because
    % fft2 and ifft2 can only work on 2D input matrices, not 3D
    % matrices, not even ones with a singleton dimension, ie., a 2D
    % matrix in disguise:
    A = squeeze(A);
    
    % Perform forward 2D FFT on GPU:
    F = fft2(A);
    
    % Process it on GPU:
    
    % First shift the spectrums origin (DC component aka 0 Hz frequency) to
    % the center of the FFT image:
    FS = fftshift(F);
    
    % Generate the amplitude spectrum, scaled to 1/100th via abs(FS)/100.0,
    % then converte the image into a texture 'fftmag' and display it:
    fftmag = GPUTypeFromToGL(1, abs(FS)/100.0, [], [], 0);
    Screen('DrawTexture', w, fftmag);
    DrawFormattedText(w, 'Amplitude spectrum post-FFT Bunny:\nPress key to continue.', 'center', 40, [0 255 0]);
    Screen('Flip', w);
    KbStrokeWait(-1);
    
    % Turn our filter 'mask' from cpu version above into a matrix on GPU:
    % We must also transpose the mask, because 'mask' was generated to
    % match the shape of the fs matrix on the cpu in Matlab/Octave. As
    % Matlab and Octave use column-major storage, whereas Psychtoolbox uses
    % row-major storage, all the GPU variables we got from Psychtoolbox
    % textures are transposed with respect to the Matlab version. A simple
    % transpose fixes this for our mask to match GPU format:
    FM = transpose(GPUsingle(mask));
    
    % Measure execution time on GPU. The cudaThreadSynchronize() command
    % makes sure we are actually measuring GPU timing with the GetSecs():
    cudaThreadSynchronize;
    tgpu = GetSecs;

    for trials = 1:10
        % Filter the amplitude spectrum by point-wise multiply with filter FM:
        FS = FM .* FS;
        
        % Shift back DC component to original position to prepare inverse FFT:
        F = ifftshift(FS);
        
        % Process inverse 2D FFT on GPU:
        B = ifft2(F);
    end
    
    cudaThreadSynchronize;
    tgpu = (GetSecs - tgpu) / trials;
    fprintf('Measured per trial runtime of GPU FFT: %f msecs [n=%i trials].\n', tgpu * 1000, trials);
    fprintf('Speedup GPU vs. CPU: %f\n', tcpu / tgpu);
    
    % Extract real component for display:
    R = real(B);
    
    % Convert R back into a floating point luminance texture 'tex' for
    % processing/display by Screen. Here it is fine to pass a
    % single-layer matrix for getting a luminance texture, because we
    % don't use Screen('MakeTexture') or similar and don't perform the
    % normalization step which would convert into a troublesome RGB
    % texture.    
    tex = GPUTypeFromToGL(1, R, [], [], 0);
    
    % Show it:
    Screen('DrawTexture', w, tex);
    DrawFormattedText(w, 'GPU-Processed FFT->Filter->IFFT Bunny:\nPress key to continue.', 'center', 40, [0 255 0]);
    Screen('Flip', w);
    KbStrokeWait(-1);
    
    % Reset interop cache before closing the windows and textures:
    GPUTypeFromToGL(4);
    
    % Close window, which will also release all textures:
    sca;
    
    % Restore sync test setting:
    Screen('Preference','SkipSyncTests', oldskip);
    
catch %#ok<CTCH>
    % Error handling.
    
    % Reset interop cache before closing the windows and textures:
    GPUTypeFromToGL(4);
    
    % Close window.
    sca;
    
    % Restore sync test setting:
    Screen('Preference','SkipSyncTests', oldskip);
    
    psychrethrow(psychlasterror);
end

% Done.
end

% Documentation of another, not yet tested, more complex and probably
% higher performance, way of doing GPU accelerated FFT with GPUmat:
%
%         fftType = cufftType;
%         fftDir  = cufftTransformDirections;
%
%         % FFT plan
%         plan = 0;
%         [status, plan] = cufftPlan2d(plan, size(d_A, 1), size(d_A, 2), fftType.CUFFT_R2C);
%         cufftCheckStatus(status, 'Error in cufftPlan2D');
%
%         % Run GPU FFT
%         [status] = cufftExecR2C(plan, getPtr(d_A), getPtr(d_B), fftDir.CUFFT_FORWARD);
%         cufftCheckStatus(status, 'Error in cufftExecR2C');
%
%         % Run GPU IFFT
%         [status] = cufftExecC2R(plan, getPtr(d_B), getPtr(d_A), fftDir.CUFFT_INVERSE);
%         cufftCheckStatus(status, 'Error in cufftExecC2C');
%
%         % results should be scaled by 1/N if compared to CPU
%         % h_B = 1/N*single(d_A);
%
%         [status] = cufftDestroy(plan);
%         cufftCheckStatus(status, 'Error in cuffDestroyPlan');
