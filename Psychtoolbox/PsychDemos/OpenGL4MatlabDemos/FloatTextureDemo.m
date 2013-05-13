function FloatTextureDemo(singlefbo, trials)

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

% Read test input image...
impath = [PsychtoolboxRoot '/PsychDemos/konijntjes1024x768.jpg']
%impath = [PsychtoolboxRoot '/PsychDemos/OpenGL4MatlabDemos/earth_512by256.jpg']
inputimage = imread(impath);
% ...convert into grayscale image...
inputimage = rgb2gray(inputimage);
% inputimage = inputimage(1:256,1:256);
% Pad it to be a square, power-of-two image:
msize = pow2(ceil(log2(max(size(inputimage)))))
tmpimage = zeros(msize, msize);
tmpimage(1:size(inputimage,1), 1:size(inputimage,2)) = inputimage(:,:);
inputimage = tmpimage;

% ...and convert it to texture:
inputtex = Screen('MakeTexture', win, inputimage);

% Retrieve an OpenGL handle to it:
texin = Screen('GetOpenGLTexture', win, inputtex);

% Setup the OpenGL rendering context of the onscreen window for use by
% OpenGL wrapper. After this command, all following OpenGL commands will
% draw into the onscreen window 'win':
Screen('BeginOpenGL', win);

% Use 
if glGetIntegerv(GL.MAX_COLOR_ATTACHMENTS_EXT)<2 || singlefbo==0
    singlefbo=0;
else
    singlefbo=1;
end;

% Create framebuffer objects and color buffer textures:
if singlefbo==0
    % Hw only supports one color attachment per FBO. Need to create two
    % FBOs for pingpong:
    [fbos(1), tex(1)]=moglCreateFBO(msize, msize);
    [fbos(2), tex(2)]=moglCreateFBO(msize, msize);
else
    % We use one FBO with two color attachments:
    [fbos(1), tex ]=moglCreateFBO(msize, msize, 2);
    % Set fbos(2) = 0 to signal not to use fbo(2).
    fbos(2)=0;
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
% instead of 0-1:
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

% Unbind FBO, reset to normal framebuffer:
moglChooseFBO(0);

% Perform reduce operation to compute texture properties:
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

% Disable our FBO's:
moglChooseFBO(0);

% Finish OpenGL rendering into PTB window and check for OpenGL errors.
Screen('EndOpenGL', win);

% Render something into visible framebuffer just to check for
% interference...
glClearColor(1,1,0,0);
glClear(GL.COLOR_BUFFER_BIT);
glColor4f(1,0,0,1);
glBegin(GL.TRIANGLES)
glVertex2d(0,0);
glVertex2d(256,0);
glVertex2d(128,256);
glEnd;
% Show rendered image at next vertical retrace:
Screen('Flip', win);

Screen('BeginOpenGL', win);

% Ok draw our image into FBO again, this time normalizing it, so that the
% range [min ; max] maps to interval [0 ; 1]:

% Draw into FBO 1:
moglChooseFBO(fbos(1), 1);

% Activate remap shader for following blit op:
glUseProgram(glslnormalizer);

minv = minv / 255;
maxv = maxv / 255;
sf = 1.0 / (maxv - minv);

% Set bias to be applied: This maps the minimum value to zero:
glUniform1f(prebias, -minv);
% No post-scale bias:
glUniform1f(postbias, 0.0);
% Scale factor to spread out everything into range 0.0 - 1.0:
glUniform1f(scalefactor, sf);

% Blit input image into it, remapping/normalizing it to range 0.0 - 1.0:
moglBlitTexture(texin);

% Disable shader:
glUseProgram(0);

% Read blitted textures back from FBOs:
% This shows the normalized 0-1 edge-map:
moglChooseFBO(fbos(1), 1);
img = glReadPixels(0, 0, msize, msize, GL_RGB, GL_FLOAT);
imshow(img);
figure

% 



if singlefbo==0
    % Dual FBO-case:
    moglChooseFBO(fbos(2), 1);
else
    % Dual color buffer on single FBO case:
    moglChooseFBO(fbos(1), 2);
end;

img2 = glReadPixels(0, 0, msize, msize, GL_RGB, GL_FLOAT);
imshow(img2, [0 255]);
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

imshow(inputimage, [0 255]);

% Well done!
return
