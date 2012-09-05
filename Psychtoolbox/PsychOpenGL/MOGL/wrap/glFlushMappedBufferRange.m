function glFlushMappedBufferRange( target, ptr, ptr )

% glFlushMappedBufferRange  Interface to OpenGL function glFlushMappedBufferRange
%
% usage:  glFlushMappedBufferRange( target, ptr, ptr )
%
% C function:  void glFlushMappedBufferRange(GLenum target, GLint ptr, GLsizei ptr)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glFlushMappedBufferRange', target, ptr, ptr );

return
