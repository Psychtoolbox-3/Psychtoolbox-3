function glInvalidateBufferSubData( buffer, ptr, ptr )

% glInvalidateBufferSubData  Interface to OpenGL function glInvalidateBufferSubData
%
% usage:  glInvalidateBufferSubData( buffer, ptr, ptr )
%
% C function:  void glInvalidateBufferSubData(GLuint buffer, GLint ptr, GLsizei ptr)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glInvalidateBufferSubData', buffer, ptr, ptr );

return
