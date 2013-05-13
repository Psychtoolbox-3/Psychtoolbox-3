function gpgpuGVF(singlefbo, trials)

if nargin < 1
    singlefbo=1;
end;

if nargin < 2
    trials = 1;
end;

% Is the script running in OpenGL Psychtoolbox?
AssertOpenGL;

% Find the screen to use for display:
screenid=max(Screen('Screens'));

% Setup Psychtoolbox for OpenGL 3D rendering support and initialize the
% mogl OpenGL for Matlab wrapper:
InitializeMatlabOpenGL(1);

% Open a double-buffered full-screen window on the main displays screen.
[win , winRect] = Screen('OpenWindow', screenid);

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Read test input image...
%impath = [PsychtoolboxRoot '/PsychDemos/konijntjes1024x768.jpg']
impath = [PsychtoolboxRoot '/PsychDemos/OpenGL4MatlabDemos/earth_512by256.jpg']
inputimage = imread(impath);
% ...convert into grayscale image...
inputimage = transpose(rgb2gray(inputimage));
%inputimage = inputimage(1:128, 1:128);
% Pad it to be a square, power-of-two image:
msize = pow2(ceil(log2(max(size(inputimage)))))
tmpimage = zeros(msize, msize);
xpad = int16(floor((msize - size(inputimage,1))/2))+1;
ypad = int16(floor((msize - size(inputimage,2))/2))+1;
tmpimage(xpad:xpad+size(inputimage,1)-1, ypad:ypad+size(inputimage,2)-1) = inputimage(:,:);
inputimage = tmpimage;

% Compute edge-map:
inputimage = uint8(edge(inputimage, 'canny'));

% TEST TEST TEST
%inputimage(:,:)=0;
%inputimage(100:105,:)=255;

Screen('EndOpenGL', win);


% ...and convert it to texture:
size(inputimage)
inputtex = Screen('MakeTexture', win, inputimage);

% Retrieve an OpenGL handle to it:
texin = Screen('GetOpenGLTexture', win, inputtex);

Screen('BeginOpenGL', win);



maxattachnr = glGetIntegerv(GL.MAX_COLOR_ATTACHMENTS_EXT)

% Use 
if glGetIntegerv(GL.MAX_COLOR_ATTACHMENTS_EXT)<2 || singlefbo==0
    singlefbo=0;
else
    singlefbo=1;
end;

precision = GL.RGBA_FLOAT16_APPLE;

% Create framebuffer objects and color buffer textures:
if singlefbo==0
    % Hw only supports one color attachment per FBO. Need to create three
    % FBOs for pingpong:
    [fbos(1), tex(1)]=moglCreateFBO(msize, msize, 1, 4, precision);
    [fbos(2), tex(2)]=moglCreateFBO(msize, msize, 1, 4, precision);
    [fbos(3), tex(3)]=moglCreateFBO(msize, msize, 1, 4, precision);
else
    % We use one FBO with two color attachments for pingpong.
    [fbos(1), tex ]=moglCreateFBO(msize, msize, 2, 4, precision);
    % Set fbos(2) = 0 to signal not to use fbo(2).
    fbos(2)=0;
    % We create one FBO for the static b, c1 and c2 arrays:
    [fbos(3), tex(3) ]=moglCreateFBO(msize, msize, 1, 4, precision);
end;

% Load our bias and rescale shader:
glslnormalizer = LoadGLSLProgramFromFiles('ScaleAndBiasShader');
prebias = glGetUniformLocation(glslnormalizer, 'prescaleoffset');
postbias = glGetUniformLocation(glslnormalizer, 'postscaleoffset');
scalefactor = glGetUniformLocation(glslnormalizer, 'scalefactor');

% Activate it for following blit op:
glUseProgram(glslnormalizer);

% Set no bias to be applied:
glUniform1f(prebias, 0.0);
glUniform1f(postbias, 0.0);

% Multiply all luminance values by 255, so they are in usual range 0-255
% instead of 0-1. We do this to reduce numeric roundoff errors.
glUniform1f(scalefactor, 255.0);

% Bind FBO 1, colorbuffer 1:
moglChooseFBO(fbos(1), 1);

% Blit input image into it, scaling it from 0-1 to 0-255:
moglBlitTexture(texin);

% Disable shader:
glUseProgram(0);

% Read it back (only the RED channel, it is a grayscale image, so
% red, green and blue contain the same content. We read it back to make
% sure we can compare Matlabs results to the GPU results:
inputimage = glReadPixels(0, 0, msize, msize, GL_RED, GL_FLOAT);
imshow(inputimage);
figure;

% Unbind FBO, reset to normal framebuffer:
moglChooseFBO(0);

% Perform reduce operation to compute image properties: We need the maximum
% and minimum luminance value:
[minv, maxv, meanv] = moglComputeMinMaxMeanOfTexture(moglGetTexForFBO(fbos(1), 1), fbos(1), fbos(2), 1);
glFinish;
tic;
for c=1:trials
    % We use the texture of fbo(1) as input image, as we know that fbo-1
    % will not be used as drawing target in first iteration. We provide
    % fbo(1) and fbo(2) as ping-pong buffers for iterative reduce
    % operation:
    [dminv, dmaxv, dmeanv] = moglComputeMinMaxMeanOfTexture(moglGetTexForFBO(fbos(1), 1), fbos(1), fbos(2), 1);
end;
durationgpu = toc / trials;

minc=min(min(inputimage));
maxc=max(max(inputimage));
meanc=mean(mean(inputimage));

tic;
minc=min(min(inputimage));
maxc=max(max(inputimage));
meanc=mean(mean(inputimage));
durationcpu = toc;

fprintf('Minimum   CPU = %f  , OpenGL = %f\n', minc, minv);
fprintf('Maximum   CPU = %f  , OpenGL = %f\n', maxc, maxv);
fprintf('Mean      CPU = %f  , OpenGL = %f\n', meanc, meanv);
fprintf('Time/pass CPU = %f ms, OpenGL = %f ms\n', durationcpu * 1000, durationgpu * 1000);


% Compute GVF in C-Code:
mu = 0.1
numiters = 100

% GVFC needs double input, not float input!
gvfcinput = double(inputimage);

tic
[gvfc_v,gvfc_u] = GVFC(gvfcinput, mu, numiters);
gvfcduration = toc * 1000
gvfcdurationperiter = gvfcduration / numiters

% This performs stage 1 of GVF initialization: Normalize all edge map
% values, remapping it so they span the whole range of 
% [min ; max] to interval [0 ; 1]. Compute initial edge map gradient by
% taking central differences as derivatives.
% We use a special shader for this that combines these two operations,
% drawing the input edge map texture into FBO 1 again:
tic;

% Load our remap, normalize and gradient map shader:
glslgradientshader = LoadGLSLProgramFromFiles('ScaleBiasAndGradientShader');
prebias = glGetUniformLocation(glslgradientshader, 'prescaleoffset');
postbias = glGetUniformLocation(glslgradientshader, 'postscaleoffset');
scalefactor = glGetUniformLocation(glslgradientshader, 'scalefactor');


% Draw into FBO 1:
moglChooseFBO(fbos(1), 1);

% Activate gradient shader for following blit op:
glUseProgram(glslgradientshader);

% Compute bias and scale to be used:
minv = minv / 255;
maxv = maxv / 255;
sf = 1.0 / (maxv - minv);

% Set bias to be applied: This maps the minimum value to zero:
glUniform1f(prebias, -minv);
% No post-scale bias:
glUniform1f(postbias, 0.0);
% Scale factor to spread out everything into range 0.0 - 1.0:
glUniform1f(scalefactor, sf);

% Blit input image into it, remapping/normalizing it to range 0.0 - 1.0 and
% then computing the gradient image on it. The RED channel will contain gx,
% the GREEN channel will contain gy, BLUE and ALPHA are unused.
moglBlitTexture(texin);

% Disable shader:
glUseProgram(0);

% Read blitted textures back from FBOs:
% This shows the normalized 0-1 edge-map:
img = glReadPixels(0, 0, msize, msize, GL_RGB, GL_FLOAT);
% Red is dx:
imshow(img(:,:,1), [-0.5 0.5]);
figure

% Green is dy:
imshow(img(:,:,2), [-0.5 0.5]);
figure

% Blue is remapped intensity:
imshow(img(:,:,3), [0 1]);
figure

% GVF Step 2: Initial flow field (u,v) is stored in (R,G) channels of
% FBO1-1's texture. Compute constant magnitude, c1 and c2 arrays in FBO 3:

% Bind FBO 3-1 for drawing. We could do fast-binding here, as dimensions
% match previously bound FBO 1-1:
moglChooseFBO(fbos(3),1);

% Load and activate our b=1-magnitude, c1, c2 shader:
glslGVFInitshader = LoadGLSLProgramFromFiles('GVFInitShader');
glUseProgram(glslGVFInitshader);

% Blit the flow field in FBO 1-1 into it:
moglBlitTexture(moglGetTexForFBO(fbos(1), 1));

% Disable shader:
glUseProgram(0);

% Read blitted textures back from FBOs:
% This shows the b=1-mag, c1 and c2 fields:
img = glReadPixels(0, 0, msize, msize, GL_RGB, GL_FLOAT);
% Red is b = 1 - mag (Range 0.5 to 1.0)
imagesc(img(:,:,1));
figure

% Green is c1:
imagesc(img(:,:,2));
figure

% Blue is c2:
imagesc(img(:,:,3));
figure

% FBO 1-1 contains initial (u,v) flow field.
% FBO 3-1 contains (b,c1,c2) static field.
% GVF iteration will bounce data between FBO 1-1 and FBO 1-2, the
% intermediate result of last iteration bound as texture to texture unit 0,
% the static field bound as texture to texture unit 1:

% Load GVF Update shader:
glslGVFUpdateshader = LoadGLSLProgramFromFiles('GVFUpdateShader');
ingradient = glGetUniformLocation(glslGVFUpdateshader, 'ingradient');
bc1c2field = glGetUniformLocation(glslGVFUpdateshader, 'bc1c2field');

% Activate update shader:
glUseProgram(glslGVFUpdateshader);

glUniform1i(ingradient, 0);
glUniform1i(bc1c2field, 1);

% Retrieve handle to mu constant and set it to reasonable value:
fourmu = glGetUniformLocation(glslGVFUpdateshader, 'fourmu');
glUniform1f(fourmu, 4 * mu);

% Bind FBO 1-2 as target for first iteration:
% We could do fast-binding as dimensions of all buffers match!
moglChooseFBO(fbos(1), 2);

% Source textures for iterations are cached in buffertex:
buffertex(1) = moglGetTexForFBO(fbos(1), 1);
buffertex(2) = moglGetTexForFBO(fbos(1), 2);

% Need to setup 2nd texture unit for constant input image:
b_c1_c2_tex = moglGetTexForFBO(fbos(3), 1);

% Unbind our constant texture:
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);

% Select the 2nd texture unit (unit 1) for setup:
glActiveTexture(GL.TEXTURE1);

% Switch it into RECTANGLE texture mapping mode:
glDisable(GL.TEXTURE_2D);
glEnable(GL.TEXTURE_RECTANGLE_EXT);

% Bind our constant texture:
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, b_c1_c2_tex);

% Make sure we use nearest neighbour sampling:
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);

% And that we clamp to edge:
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);

% Choose texture application function to be a neutral REPLACE:
glTexEnvfv(GL.TEXTURE_ENV,GL.TEXTURE_ENV_MODE,GL.REPLACE);

% Unit 1 is ready, switch back to Unit 0, the first unit:
glActiveTexture(GL.TEXTURE0);

% Switch it into RECTANGLE texture mapping mode:
glDisable(GL.TEXTURE_2D);
glEnable(GL.TEXTURE_RECTANGLE_EXT);

inputtex = buffertex(1);

glFinish;

gvfinitduration = toc * 1000

tic;

% GVF iterative ping-pong update loop:
bufferid = 0;
for i=1:numiters
%tic
    % Compute drawbuffers id and its texturehandle:
    bufferid = 1 - bufferid;
    futuretex = buffertex(bufferid + 1);

    % Bind proper FBO as target for intermediate reduce results:
    if singlefbo == 0
        % Use two FBO's for pingpong:
        if bufferid==0
            glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, fbos(1));
        else
            glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, fbos(2));
        end;
    else
        % Use one dual-buffer FBO for pingpong:
        % We only switch the draw-buffer to keep the overhead as low as
        % possible:
        glDrawBuffer(GL.COLOR_ATTACHMENT0_EXT + bufferid);
    end;
    
    % Do it: We use fast-blit mode (1) because we did texture setup
    % ourselves:
    moglBlitTexture(inputtex, 0, 0, msize, msize, 1);

    % Assign source texture for next pass:
    inputtex = futuretex;
%glFinish;
%gpupass = toc * 1000
end;

glFinish
gvfduration = toc * 1000
gvfperiteration = gvfduration / numiters

% Shut down the GVF iteration shader...
glUseProgram(0);
% ...and both texture units:
glDisable(GL.TEXTURE_RECTANGLE_EXT);
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
glActiveTexture(GL.TEXTURE1);
glDisable(GL.TEXTURE_RECTANGLE_EXT);
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
glActiveTexture(GL.TEXTURE0);

% Assign read buffer for GVF readout:
glReadBuffer(GL.COLOR_ATTACHMENT0_EXT + bufferid);

% Readout final result: RED = flow_x, GREEN = flow_y:
gvf_finalimg = glReadPixels(0, 0, msize, msize, GL_RGB, GL_FLOAT);

% Final flow_x:
imagesc(gvf_finalimg(:,:,1));
figure

% Final flow_y:
imagesc(gvf_finalimg(:,:,2));
figure

% Unbind FBO, reset to normal framebuffer:
moglChooseFBO(0);

% Delete all FBO's:
for i=1:length(fbos)
    moglDeleteFBO(fbos(i));
end;

% Shut down OpenGL rendering:
Screen('EndOpenGL', win);

% Close onscreen window and release all other ressources:
Screen('CloseAll');

%minu = min(min(gvfc_u))
minu = min(min(gvf_finalimg(:,:,1)))
%maxu = max(max(gvfc_u))
maxu = max(max(gvf_finalimg(:,:,1)))

return

imagesc(gvfc_u);
figure;
imagesc(gvfc_v);
figure;
diff_u = gvfc_u - gvf_finalimg(:,:,1);
diff_v = gvfc_v - gvf_finalimg(:,:,2);
imagesc(diff_u);
figure;
imagesc(diff_v);

maxdiff=max(max(abs(diff_u) + abs(diff_v)))
avgdiff=mean(mean(abs(diff_u) + abs(diff_v)))
% Well done!
return
