function glCopyBufferSubData( readTarget, writeTarget, ptr, ptr, ptr )

% glCopyBufferSubData  Interface to OpenGL function glCopyBufferSubData
%
% usage:  glCopyBufferSubData( readTarget, writeTarget, ptr, ptr, ptr )
%
% C function:  void glCopyBufferSubData(GLenum readTarget, GLenum writeTarget, GLint ptr, GLint ptr, GLsizei ptr)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glCopyBufferSubData', readTarget, writeTarget, ptr, ptr, ptr );

return
