function CreateResolutionPyramid(tex, glsl, usebilinear)
% Build a mip-map image resolution pyramid for given texture.
%
% CreateResolutionPyramid(tex [, glsl=0][, usebilinear=0])
%
% 'tex' must be the Screen() texture or offscreen window handle of a
% texture or offscreen window of type GL_TEXTURE_2D, ie., created with
% 'specialFlags' setting 1. The function will create a OpenGL mip-map
% resolution pyramid by successive downsampling of 'tex' image content down
% to a image of 1x1 pixels. If 'glsl' is provided as handle to a shader, it
% will use that shader to do the downsampling, otherwise it will use
% bilinear filtering for successive pyramid levels. If 'usebilinear' is set
% to 1 it will feed the shader with linear interpolated samples, otherwise
% with nearest neighbour samples.
%
% You will likely also want to set 'specialFlags' flag 8 or 16 to prevent
% Screen('DrawTexture') from overwriting the image pyramid with its own
% auto-generated one.
%
% This needs a modern GLSL shading capable graphics card to work.
%

% History:
% 27.08.2012  mk   Written.
%

% Needed for GL constant definitions:
global GL;

% Make sure MOGL is initialized:
if isempty(GL)
    InitializeMatlabOpenGL([],[],1);
end

if nargin < 1 || isempty(tex)
    error('Must provide a valid Screen() texture or offscreen window handle for building of resolution pyramid!');
end

if nargin < 2 || isempty(glsl)
    glsl = 0;
end

if nargin < 3 || isempty(usebilinear)
    usebilinear = 0;
end

% Only supported for offscreen targets, not onscreen windows, proxies etc.:
if Screen('WindowKind', tex) ~= -1
    error('Tried to built a resolution image pyramid for something else than a texture or offscreen window! Unsupported.');
end

% This, as a well defined side effect, will bind the tex'tures associated
% FBO, or create one if it doesn't have one yet, as well as normalize the
% textures format and orientation, just as we need it:
Screen('GetWindowInfo', tex);

% Get width and height:
[w, h] = Screen('WindowSize', tex);

% Get low-level texture handle and type:
[gltexid, gltextarget] = Screen('GetOpenGLTexture', tex, tex);

% Make sure it is a mip-map capable type:
if gltextarget ~= GL.TEXTURE_2D
    error('Tried to built a resolution image pyramid for something else than a GL_TEXTURE_2D texture! Unsupported.');
end

% Time for a state backup:
glPushAttrib(GL.ALL_ATTRIB_BITS);

% Ok, the currently bound FBO is the FBO to use for drawing into this
% texture. And we got the gltexid of the corresponding texture.
% We also know that all texture parameters are set in a way to make it
% available as a rendertarget, e.g., the texture minification/magnification
% filters, wrap modes etc. are set to safe values, e.g., nearest neighbour
% sampling, clamp to edge etc.

% Bind the texture:
glBindTexture(gltextarget, gltexid);

% Enable sampling:
glDisable(GL.TEXTURE_RECTANGLE_EXT);
glEnable(gltextarget);

% Compute to what miplevel we need to downsample to cover all bases:
maxlod = floor(log2(max(w, h)));

% Query size of mip-level 'maxlod':
tw = glGetTexLevelParameteriv(gltextarget, maxlod, GL.TEXTURE_WIDTH);
th = glGetTexLevelParameteriv(gltextarget, maxlod, GL.TEXTURE_HEIGHT);

% If no such mip-level exists?
if tw == 0 || th == 0
    % No such mip-level: Generate initial mip-chain for this texture. We
    % allocate mip-levels of proper format and size, but empty, ie.,
    % without actually computing any content here:
    internalFormat = glGetTexLevelParameteriv(gltextarget, 0, GL.TEXTURE_INTERNAL_FORMAT);

    for miplevel = 1:maxlod
        % Need to case tw, th to double(), use max() operators and whatnot
        % to get proper size specs -- Matlab arithmetic on the returned
        % int32's from query behaves rather strange otherwise:
        tw = double(glGetTexLevelParameteriv(gltextarget, miplevel - 1, GL.TEXTURE_WIDTH));
        th = double(glGetTexLevelParameteriv(gltextarget, miplevel - 1, GL.TEXTURE_HEIGHT));
        glTexImage2D(gltextarget, miplevel, internalFormat, max(floor(tw/2), 1), max(floor(th/2), 1), 0, GL.RGBA, GL.UNSIGNED_BYTE, 0);
    end
end

% Bind downsampling shader and set it up to read from texture unit 0:
glUseProgram(glsl);

% If we've actually bound the fixed-function pipeline, we skip this step:
if glsl ~= 0
    glUniform1i(glGetUniformLocation(glsl, 'Image'), 0);
    
    % Query locations of uniforms for passing src / dst dimensions:
    srcSize = glGetUniformLocation(glsl, 'srcSize');
    dstSize = glGetUniformLocation(glsl, 'dstSize');
end

% Just to be extra safe:
glTexParameteri(gltextarget, GL.TEXTURE_WRAP_S, GL.CLAMP_TO_EDGE);
glTexParameteri(gltextarget, GL.TEXTURE_WRAP_T, GL.CLAMP_TO_EDGE);

if glsl ~= 0 && ~usebilinear
    % Use nearest neighbour sampling for real downsampling shaders, so they
    % have full control over actual sampling:
    glTexParameteri(gltextarget, GL.TEXTURE_MIN_FILTER, GL.NEAREST);
    glTexParameteri(gltextarget, GL.TEXTURE_MAG_FILTER, GL.NEAREST);
else
    % For fixed function, use bilinear filtering. This should pretty much
    % emulate a 2x2 texels box-filter:
    glTexParameteri(gltextarget, GL.TEXTURE_MIN_FILTER, GL.LINEAR);
    glTexParameteri(gltextarget, GL.TEXTURE_MAG_FILTER, GL.LINEAR);
end

% Base level 0 is always enabled, as we never write to it:
glTexParameteri(gltextarget, GL.TEXTURE_BASE_LEVEL, 0);

% Iterate over all mip-levels 'miplevel' from 1 to maxlod, and source each
% 'miplevel-1' to compute 'miplevel':
for miplevel = 1:maxlod
    % Query size of target mip level:
    tw = glGetTexLevelParameteriv(gltextarget, miplevel, GL.TEXTURE_WIDTH);
    th = glGetTexLevelParameteriv(gltextarget, miplevel, GL.TEXTURE_HEIGHT);

    % If size is smaller than 1 texel in any dimension, exit the downsampling loop:
    if tw < 1 || th < 1
        break;
    end
    
    % Set mip-level upper range limit for sampling from texture from to "miplevel-1":
    glTexParameteri(gltextarget, GL.TEXTURE_MAX_LEVEL, miplevel-1);
    
    % Select the target mip-level for drawing into (== render-to-texture
    % via FBO color attachment) to 'miplevel':
    glFramebufferTexture2DEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, gltextarget, gltexid, miplevel);

    % Setup drawing transforms, projection and viewport:
    glMatrixMode(GL.PROJECTION);
    glLoadIdentity;
    gluOrtho2D(0, tw, 0 , th);
    glMatrixMode(GL.MODELVIEW);
    glLoadIdentity;
    glViewport(0, 0, tw, th);
    
    % Ok, we should be setup properly now. Any shader sampling from the
    % texture can only do nearest neighbour sampling from level
    % 'miplevel-1', and all writes go to 'miplevel'.
    
    % Communicate sizes of source and destination surface to possible
    % shaders:
    if glsl
        glUniform2f(srcSize, glGetTexLevelParameteriv(gltextarget, miplevel-1, GL.TEXTURE_WIDTH), glGetTexLevelParameteriv(gltextarget, miplevel-1, GL.TEXTURE_HEIGHT));
        glUniform2f(dstSize, tw, th);
    end
    
    % Blit a "fullscreen quad" - or more specifically a full-texture quad
    % to drive the downsampling shader:
    glColor4f(1,1,1,1);
    glBegin(GL.QUADS)
    glTexCoord2f(0,0);
    glVertex2f(0,0);
    glTexCoord2f(0,1);
    glVertex2f(0,th);
    glTexCoord2f(1,1);
    glVertex2f(tw,th);
    glTexCoord2f(1,0);
    glVertex2f(tw,0);
    glEnd;
    % Rinse, wash, repeat at next miplevel...
end

% Reset mip-level range to normal:
glTexParameteri(gltextarget, GL.TEXTURE_BASE_LEVEL, 0);
glTexParameteri(gltextarget, GL.TEXTURE_MAX_LEVEL, 1000);
    
% Disable and unbind texture:
glDisable(gltextarget);
glBindTexture(gltextarget, 0);

% Restore framebuffer attachment for texture to render into miplevel zero:
glFramebufferTexture2DEXT(GL.FRAMEBUFFER_EXT, GL.COLOR_ATTACHMENT0_EXT, gltextarget, gltexid, 0);

% Unbind shader:
glUseProgram(0);

% Final state: Texture unbound from sampler, but bound as current
% framebuffer/rendertarget. Regular PTB code will take care of this.

% Restore pre-op state:
glPopAttrib;

% Done.
return;
