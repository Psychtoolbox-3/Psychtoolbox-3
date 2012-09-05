function img = glGetnCompressedTexImageARB( target, lod, bufSize )

% glGetnCompressedTexImageARB  Interface to OpenGL function glGetnCompressedTexImageARB
%
% usage:  img = glGetnCompressedTexImageARB( target, lod, bufSize )
%
% C function:  void glGetnCompressedTexImageARB(GLenum target, GLint lod, GLsizei bufSize, GLvoid* img)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

% Query size of compressed image and alloc properly sized output buffer:
img = uint8(zeros(1, bufSize));

moglcore( 'glGetnCompressedTexImageARB', target, lod, bufSize, img );

return
