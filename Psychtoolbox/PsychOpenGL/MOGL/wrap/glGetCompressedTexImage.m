function img = glGetCompressedTexImage( target, lod )

% glGetCompressedTexImage  Interface to OpenGL function glGetCompressedTexImage
%
% usage:  img = glGetCompressedTexImage( target, lod )
%
% C function:  void glGetCompressedTexImage(GLenum target, GLint lod, GLvoid* img)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

global GL;

if nargin~=2,
    error('invalid number of arguments');
end

% Query size of compressed image and alloc properly sized output buffer:
img = uint8(zeros(1, glGetTexLevelParameteriv(target, lod, GL.TEXTURE_COMPRESSED_IMAGE_SIZE)));

% Go get it:
moglcore( 'glGetCompressedTexImage', target, lod, img );

return
