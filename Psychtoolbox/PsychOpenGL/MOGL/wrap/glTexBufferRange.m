function glTexBufferRange( target, internalformat, buffer, ptr, ptr )

% glTexBufferRange  Interface to OpenGL function glTexBufferRange
%
% usage:  glTexBufferRange( target, internalformat, buffer, ptr, ptr )
%
% C function:  void glTexBufferRange(GLenum target, GLenum internalformat, GLuint buffer, GLint ptr, GLsizei ptr)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glTexBufferRange', target, internalformat, buffer, ptr, ptr );

return
