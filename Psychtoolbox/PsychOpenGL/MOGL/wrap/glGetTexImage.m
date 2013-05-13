function pixels = glGetTexImage( target, level, format, type, bufSize )

% glGetTexImage  Interface to OpenGL function glGetTexImage
%
% usage:  pixels = glGetTexImage( target, level, format, type )
%
% C function:  void glGetTexImage(GLenum target, GLint level, GLenum format, GLenum type, GLvoid* pixels)

% 05-Mar-2006 -- created (generated automatically from header files)
% 05-Sep-2006 -- Made functional by implementing buffer allocation code. (MK)

% ---allocate---
% ---protected---
% ---skip---

if nargin~=4 && nargin~=5,
    error('invalid number of arguments');
end

% check format and type
global GL

switch(format)
    case {GL.RED, GL.GREEN, GL.BLUE, GL.ALPHA, GL.LUMINANCE, GL.INTENSITY, GL.DEPTH_COMPONENT, GL.STENCIL_INDEX }
        numperpixel = 1;
    case {GL.LUMINANCE_ALPHA }
        numperpixel = 2;
    case {GL.RGB , GL.BGR }
        numperpixel = 3;
    case {GL.RGBA , GL.BGRA }
        numperpixel = 4;
    otherwise
        error('Invalid format passed to glGetTexImage.');
end;    
        
switch(type)
    case {GL.UNSIGNED_BYTE , GL.BYTE }
        numsize = 1;
    case { GL.UNSIGNED_SHORT , GL.SHORT }
        numsize = 2;
    case { GL.UNSIGNED_INT , GL.INT , GL.FLOAT }
        numsize = 4;
    otherwise
        error('Invalid type passed to glGetTexImage.');
end;

% Query size of currently bound texture object:
width = glGetTexLevelParameteriv(target, level, GL.TEXTURE_WIDTH);
height = glGetTexLevelParameteriv(target, level, GL.TEXTURE_HEIGHT);

if width <=0 || height<=0
    error('Invalid (negative or zero) size of texture image object.');
end

% Allocate memory:
pixels=zeros(numperpixel, width, height);

% Tell OpenGL that we accept byte-aligned aka unaligned data.
glPixelStorei(GL.PACK_ALIGNMENT, 1);

% Perform proper type-cast:
switch(type)
    case GL.UNSIGNED_BYTE
        pixels = uint8(pixels);
    case GL.BYTE
        pixels = int8(pixels);
    case GL.UNSIGNED_SHORT
        pixels = uint16(pixels);
    case GL.SHORT
        pixels = int16(pixels);
    case GL.UNSIGNED_INT
        pixels = uint32(pixels);
    case GL.INT
        pixels = int32(pixels);
    case GL.FLOAT
        pixels = single(pixels);
end;

if exist('bufSize', 'var')
    % This is actually a call to glGetnTexImageARB:
    moglcore( 'glGetnTexImageARB', target, level, format, type, bufSize, pixels );
else
    % Execute actual call to moglcore:
    moglcore( 'glGetTexImage', target, level, format, type, pixels );
end

return
