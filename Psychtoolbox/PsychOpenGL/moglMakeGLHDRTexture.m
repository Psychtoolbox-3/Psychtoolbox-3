function gltexId = moglMakeGLHDRTexture(hdrImage, gltextarget, halffloat)
% gltexId = moglMakeGLHDRTexture(hdrImage, gltextarget [, halffloat])
%
% Create a high dynamic range OpenGL texture from the given
% Matlab/Octave high dynamic range image matrix 'hdrImage' and attach
% it to the OpenGL texture target of type 'gltextarget'.
%
% Before calling this function you must have called Screen('BeginOpenGL', win);
% with 'win' being the window handle for the onscreen window to
% which the texture should be attached.
%
% Returns an OpenGL texture handle 'gltexId' for the created texture.
% The texture can be used like any other OpenGL texture, just
% with the difference that it represents its color values with 32 bit
% floating point precision instead of 8 bpc fixed point precision.
%
% 'hdrImage' must be a (height, width, 4) matrix of type 'double' or
% 'single', where channel 1=Red, 2=Green, 3=Blue, 4=Alpha. The numeric
% range (0.0 - 1.0) maps to (minimum intensity - maximum intensity).
%
% If 'halffloat' is set to 1 (default is 0), then the texture is created as
% half-precision floating point texture, ie 16 bpc floats. This saves
% memory and bandwidth and allows for floating point filtering on Geforce
% 6000 and 7000 series.
%
% If you want to create a pure 2D texture and use it with Psychtoolbox'
% standard Screen() drawing command, then call moglMakeHDRTexture() instead.
%
% This function is most useful if you want to create textures not handled
% by Psychtoolbox, i.e., cube map textures for environment mapping and
% texture based lighting.

% History:
% 7.8.2006 Created (MK)

persistent alreadychecked;
global GL;

if isempty(GL)
   InitializeMatlabOpenGL;
end;

if isempty(alreadychecked)
   % Check if required OpenGL extensions are available.
   try
      extensions = glGetString(GL.EXTENSIONS);
   catch
      error('MakeGLHDRTexture called before caling Screen(''BeginOpenGL'') on an onscreen window! This will not work...');
   end;

   if isempty(findstr(extensions, 'GL_ARB_texture_float'))
      error('This combination of driver and graphics hardware does not support floating point textures as required for HDR.');
   end;

   alreadychecked = 1;
end

if nargin < 1
   error('Required input ''hdrImage'' matrix is missing.');
end

if nargin < 2
   error('Required OpenGL texture target ''gltextarget'' is missing.');
end

if nargin < 3
    halffloat = 0;
end

if size(hdrImage,1)<1 || size(hdrImage,2)<1
   error('Passed ''hdrImage'' does not have required minimum width x height of at least 1 by 1 pixels.');
end

if size(hdrImage,3)~=4
   error('Passed ''hdrImage'' does not have required number of 4 layers (R,G,B,A).');
end

% Convert into single precision floating-point, whatever it was before.
hdrImage = moglsingle(hdrImage);

% Transpose it and merge from multi-layer to multi-component-per-pixel format:
% Also flip it vertically to match OpenGL style.
for i=1:4
  glImage(i,:,:)=transpose(flipud(hdrImage(:,:,i)));
end

% Build a proper OpenGL texture:

% Generate texture object:
gltexId=glGenTextures(1);

% Hack, needs to be improved...
if IsLinux
    % Redefine vor NVidia:
    GL.RGBA_FLOAT32_APPLE = hex2dec('8883');
end;

% Bind the new texture object:
glBindTexture(gltextarget, gltexId);

if halffloat
    intformat = GL.RGBA_FLOAT16_APPLE;
else
    intformat = GL.RGBA_FLOAT32_APPLE;
end

% Assign the transformed HDR image as mipmap level zero:
glTexImage2D(gltextarget, 0, intformat, size(hdrImage,2), size(hdrImage,1), 0, GL.RGBA, GL.FLOAT, glImage);

% Setup texture wrapping behaviour to clamp, as other behaviours are
% unsupported on many gfx-cards for rectangle textures:
glTexParameterfv(gltextarget,GL.TEXTURE_WRAP_S,GL.CLAMP);
glTexParameterfv(gltextarget,GL.TEXTURE_WRAP_T,GL.CLAMP);

% Setup filtering for the texture - Use nearest neighbour as floating
% point filtering usually unsupported.
glTexParameterfv(gltextarget,GL.TEXTURE_MAG_FILTER,GL.NEAREST);
glTexParameterfv(gltextarget,GL.TEXTURE_MIN_FILTER,GL.NEAREST);

% Unbind texture after setup:
glBindTexture(gltextarget, 0);

% Texture should be ready for use now. Return its handle in 'gltexId':
return;
