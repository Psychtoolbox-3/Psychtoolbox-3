function glDebugMessageInsertARB( source, type, id, severity, length, buf )

% glDebugMessageInsertARB  Interface to OpenGL function glDebugMessageInsertARB
%
% usage:  glDebugMessageInsertARB( source, type, id, severity, length, buf )
%
% C function:  void glDebugMessageInsertARB(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* buf)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glDebugMessageInsertARB', source, type, id, severity, length, uint8(buf) );

return
