function retpixels = glReadnPixelsARB( x, y, width, height, format, type, bufSize )

% glReadnPixelsARB  Interface to OpenGL function glReadnPixelsARB
%
% usage:  data = glReadnPixelsARB( x, y, width, height, format, type, bufSize )
%
% C function:  void glReadnPixelsARB(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLsizei bufSize, GLvoid* data)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

global GL;

if nargin~=7,
    error('invalid number of arguments');
end

% check format and type

if x < 0 || y < 0
    error('Invalid (negative) (x,y) offset passed to glReadPixels.');
end

if width <= 0 || height<=0
    error('Invalid (negative or zero) (width, height) passed to glReadPixels.');
end

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
        error('Invalid format passed to glReadPixels.');
end

% Readback to host memory, aka Matlab- or Octave- matrix:
% Allocate memory:
pixels=zeros(numperpixel, width, height);

% Tell OpenGL that we accept byte-aligned aka unaligned data.
glPixelStorei(GL.PACK_ALIGNMENT, 1);

% Perform proper type-cast:
switch(type)
    case GL.UNSIGNED_BYTE
        pixels = uint8(pixels);
        pclass = 'uint8';
    case GL.BYTE
        pixels = int8(pixels);
        pclass = 'int8';
    case GL.UNSIGNED_SHORT
        pixels = uint16(pixels);
        pclass = 'uint16';
    case GL.SHORT
        pixels = int16(pixels);
        pclass = 'int16';
    case GL.UNSIGNED_INT
        pixels = uint32(pixels);
        pclass = 'uint32';
    case GL.INT
        pixels = int32(pixels);
        pclass = 'int32';
    case GL.FLOAT
        pixels = moglsingle(pixels);
        pclass = 'double';
    otherwise
        error('Invalid type argument passed to glReadPixels()!');
end;

% Execute actual call:
moglcore( 'glReadnPixelsARB', x, y, width, height, format, type, bufSize, pixels );

% Rearrange data in Matlab friendly format:
retpixels = zeros(size(pixels,2), size(pixels,3), size(pixels,1), pclass);
for i=1:numperpixel
    retpixels(:,:,i) = pixels(i,:,:);
end
return
