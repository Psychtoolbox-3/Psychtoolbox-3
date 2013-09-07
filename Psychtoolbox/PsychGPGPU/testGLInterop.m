function testGLInterop

PsychDefaultSetup(1);

GPUstart;
%memcpyCudaOpenGL(6, 1); % Be verbose.
%memcpyCudaOpenGL(5, 50); % Use 50 slots in LRU cache.

Screen('Preference','SkipSyncTests', 2);
PsychImaging('PrepareConfiguration');
%PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
w = PsychImaging('OpenWindow', 0, 0);

tex = ones(200, 200, 4);

% Read our beloved bunny image from filesystem:
bunnyimg = imread([PsychtoolboxRoot 'PsychDemos/konijntjes1024x768.jpg']);
bunnyimg(:,:,4) = 255;

% Maketexture in float precision with upright texture orientation:
bunnytex = Screen('MakeTexture', w, double(bunnyimg(:,:,1:2))/255, [], [], 2, 1);

bunre = Screen('Rect', bunnytex)
bunp  = Screen('Pixelsize', bunnytex)

%bunnytex = Screen('Openoffscreenwindow', w, [1 1 0 0], [0 0 1024 768], 128);
%DrawFormattedText(bunnytex, 'HELLO WORLD!', 'center', 'center', [255 0 0]);

t = bunnytex;

Screen('DrawTexture', w, t);

if 0
    T = GPUTypeFromToGL(0, w, 2);
end
Screen('Flip', w);


if 1
    if 1
        % Test high-level interface for pure mortals:
        T = GPUTypeFromToGL(0, bunnytex, [], [], 0);
    else
        % Test low-level OpenGL object interface which does without any calls
        % into Screen() -- Important for interop with non-ptb code and for use
        % within imaging pipeline, ie., within Screen's execution context, as
        % Screen is currently non-reentrant.
        [texstruct.glhandle, texstruct.gltarget] = Screen('GetOpenGLTexture', w, t);
        [texstruct.width, texstruct.height] = Screen('Windowsize', t);
        texstruct.bpp = Screen('Pixelsize', t) / 8;
        texstruct.nrchannels = 4;
        T = GPUTypeFromToGL(0, texstruct, 1);
    end
end

rect = Screen('Rect', t)

%T = GPUsingle(rand(1, 1024, 768));
foo = size(T)
H = T;
foo = size(H)
c = 0;
t0 = GetSecs;
while c < 500
    %T = T .* 0.99;
    H = T .* (0.5 + 0.5 * sin(GetSecs * 10));
    if 1
        t = GPUTypeFromToGL(1, H, 0, t, 0);
        %ptr = GPUTypeFromToGL(6, t, 0);
            %GPUTypeFromToGL(3, t, 0);
            % Screen('Close', t);
            % t = GPUTypeFromToGL(1, H, 0, [], 1);
        
        Screen('DrawTexture', w, t, [], [], [], 0);
    else
        GPUTypeFromToGL(1, H, 2, w);
    end
    Screen('Flip', w, 0, 2, 2);
    c = c + 1;
end
Screen('DrawingFinished', w, 2, 1);
fps = c / (GetSecs - t0);
gbs = fps * 2 * numel(T) * 4 / 1024 / 1024 / 1024;
fprintf('FPS = %f [%f GB/sec]\n', fps, gbs);
GPUTypeFromToGL(4);
sca;

