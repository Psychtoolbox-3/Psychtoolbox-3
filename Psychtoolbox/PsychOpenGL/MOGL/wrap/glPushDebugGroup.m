function glPushDebugGroup( source, id, length, message )

% glPushDebugGroup  Interface to OpenGL function glPushDebugGroup
%
% usage:  glPushDebugGroup( source, id, length, message )
%
% C function:  void glPushDebugGroup(GLenum source, GLuint id, GLsizei length, const GLchar* message)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glPushDebugGroup', source, id, length, uint8(message) );

return
