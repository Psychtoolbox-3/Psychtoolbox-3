function glClearBufferSubData( target, internalformat, ptr, ptr, format, type, data )

% glClearBufferSubData  Interface to OpenGL function glClearBufferSubData
%
% usage:  glClearBufferSubData( target, internalformat, ptr, ptr, format, type, data )
%
% C function:  void glClearBufferSubData(GLenum target, GLenum internalformat, GLint ptr, GLsizei ptr, GLenum format, GLenum type, const void* data)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glClearBufferSubData', target, internalformat, ptr, ptr, format, type, data );

return
