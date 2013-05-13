function [minv, maxv, meanv] = moglComputeMinMaxMeanOfTexture(texid, pingpongfbo1, pingpongfbo2, fastpath)
% [minv, maxv, meanv] = moglComputeMinMaxMeanOfTexture(texid, pingpongfbo1, pingpongfbo2, fastpath)
%
% This function expects a square RGB texture as input, whose color channels
% contain identical values L=R=G=B where L is a grayscale (luminance) input
% image. It computes and returns its global minimum, maximum and mean luminance
% on the GPU.
%
% The width and height of the texture need to be divisible by two (even numbers).
% It only works on RECTANGLE textures, not on standard power-of-two textures!
% You are not allowed to use an input texture 'texid' that is identical to
% the color buffer texture of the FBO 'pingpongfbo2'!

% Notes:
% FBO ping-ponging takes 13.7 ms on GeForce-7800 for 256 by 256 image.
%
% History:
% 30.5.2006 Written (MK).

% Need handle to GL constant definitions:
global GL;

% We cache our reduction shader for reuse on later invocations:
persistent reduceshader;
if isempty(reduceshader)
    reduceshader = LoadGLSLProgramFromFiles('MinMaxMeanReduceShader');
end;

if nargin < 4
    fastpath = 0;
else
    if isempty(fastpath)
        fastpath = 0;
    end;
end;

if nargin < 2
    error('Need texid of input texture and the fbo-handles of at least the first ping-pong FBO!');
end;

% Query size of input texture:
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, texid);
w = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_WIDTH);
h = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_HEIGHT);

% We skip all error checking and such if fastpath > 0. This saves about 4
% ms of setup time.
if fastpath == 0
    % Take the slow, but safe route:    
    if ~glIsFramebufferEXT(pingpongfbo1)
        error('Invalid pingpongfbo1 identifier passed. This is not a valid FBO.');
    end;

    if nargin < 3
        pingpongfbo2 = 0;
    end;

    if isempty(pingpongfbo2)
        pingpongfbo2 = 0;
    end;

    if pingpongfbo2 > 0
        % Second bounce buffer provided. Check it.
        if ~glIsFramebufferEXT(pingpongfbo2)
            error('Invalid pingpongfbo2 identifier passed. This is not a valid FBO.');
        end;
        singlefbo = 0;
    else
        % Only one FBO provided, make sure that one has two color attachments:
        glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo1);
        if glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT1_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE_EXT)~=GL.TEXTURE
            error('Only pingpongfbo1 provided, not pingpongfbo2 *and* pingpongfbo1 does not have 2 color buffers. This will not work...');
        end;
        singlefbo = 1;
    end;

    if w~=h || w<2 || mod(w,2)>0
        error('Invalid size of input texture: Needs to be even, width = height and at least 2 pixels wide.');
    end;

    % Query size of FBO's:
    glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo1);
    fbotexid = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, fbotexid);
    fw = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_WIDTH);
    fh = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_HEIGHT);
    glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
    if fw < w/2 || fh < h/2
        error('Pingpong FBO 1 is too small for reduce operation!');
    end;

    % Query size of FBO's:
    if singlefbo == 0
        glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo2);
        fbotexid = glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
        glBindTexture(GL.TEXTURE_RECTANGLE_EXT, fbotexid);
        fw = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_WIDTH);
        fh = glGetTexLevelParameteriv(GL.TEXTURE_RECTANGLE_EXT, 0, GL.TEXTURE_HEIGHT);
        glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
        if fw < w/2 || fh < h/2
            error('Pingpong FBO 2 is too small for reduce operation!');
        end;
    end;

    % Make sure we don't try to read from the pingpongfbo2's texture as input
    % texture - that would cause undefined behaviour.
    if singlefbo == 0
        if texid == moglGetTexForFBO(pingpongfbo2)
            error('Input texture texid is the texture attached to pingpongfbo2! Thats forbidden! Aborted.');
        end;
    else
        if texid == moglGetTexForFBO(pingpongfbo1, 2)
            error('Input texture texid is the texture attached to color buffer 2 of pingpongfbo1! That is forbidden! Aborted.');
        end;
    end;
else
    % Fast path: Minimal setup, no error checks:
    if pingpongfbo2 > 0
        singlefbo = 0;
    else
        singlefbo = 1;
    end;
end;

% Ok we survived error-checking, everything is fine. Do it:
reducestep = 2;
numiters = ceil(log(double(w))/log(reducestep));
currentsize = w;

% On first iteration we operate on the input texture:
inputtex = texid;

% Setup texture mapping:
glDisable(GL.TEXTURE_2D);
glEnable(GL.TEXTURE_RECTANGLE_EXT);

% Make sure we use nearest neighbour sampling:
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
% And that we clamp to edge:
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_S, GL.CLAMP);
glTexParameteri(GL.TEXTURE_RECTANGLE_EXT, GL.TEXTURE_WRAP_T, GL.CLAMP);

% Do a full setup of the initial FBO and colorbuffer before iterating,
% because this is very costly. The iteration loops only do the minimal
% amount of switching via glXXX calls to reduce switch-overhead:
if singlefbo==1
    % One FBO, two buffers:
    moglChooseFBO(pingpongfbo1, 2);
    buffertex(1)=glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
    buffertex(2)=glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT1_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);
else
    % Two FBO's:
    glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo1);
    buffertex(1)=glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);    
    moglChooseFBO(pingpongfbo2, 1);
    buffertex(2)=glGetFramebufferAttachmentParameterivEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, GL.FRAMEBUFFER_ATTACHMENT_OBJECT_NAME_EXT);    
end;

% Activate reduce-shader:
glUseProgram(reduceshader);

% Do numiters pingpong-reduce operations:
bufferid = 0;
for i=1:numiters
    % Compute drawbuffers id and its texturehandle:
    bufferid = 1 - bufferid;
    futuretex = buffertex(bufferid + 1);

    % Bind proper FBO as target for intermediate reduce results:
    if singlefbo == 0
        % Use two FBO's for pingpong:
        if bufferid==0
            glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo1);
        else
            glBindFramebufferEXT(GL.FRAMEBUFFER_EXT, pingpongfbo2);
        end;
    else
        % Use one dual-buffer FBO for pingpong:
        % We only switch the draw-buffer to keep the overhead as low as
        % possible:
        glDrawBuffer(GL.COLOR_ATTACHMENT0_EXT + bufferid);
    end;
    
    % Perform reduction blit operation to a size of 'currentsize' by 'currentsize':
    currentsize = floor(currentsize / reducestep);

    % Do it: We use fast-blit mode (1) because we did texture setup
    % ourselves:
    moglBlitTexture(inputtex, 0, 0, currentsize, currentsize, 1);

    % Assign source texture for next pass:
    inputtex = futuretex;
end;

% Ok, pixel (0,0) of our currently bound FBO should contain the final
% result. Read it out:
if singlefbo == 0, bufferid = 0; end;
glReadBuffer(GL.COLOR_ATTACHMENT0_EXT + bufferid);
resultpixel = double(glReadPixels(0, 0, 1, 1, GL.RGB, GL.FLOAT));
meanv = resultpixel(1);
maxv  = resultpixel(2);
minv  = resultpixel(3);

% Disable texture mapping...
glDisable(GL.TEXTURE_RECTANGLE_EXT);
glBindTexture(GL.TEXTURE_RECTANGLE_EXT, 0);
% ...disable reduce - shader:
glUseProgram(0);

% Well done!
return;
