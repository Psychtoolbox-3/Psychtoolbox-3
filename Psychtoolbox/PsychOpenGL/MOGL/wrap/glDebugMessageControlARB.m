function glDebugMessageControlARB( source, type, severity, count, ids, enabled )

% glDebugMessageControlARB  Interface to OpenGL function glDebugMessageControlARB
%
% usage:  glDebugMessageControlARB( source, type, severity, count, ids, enabled )
%
% C function:  void glDebugMessageControlARB(GLenum source, GLenum type, GLenum severity, GLsizei count, const GLuint* ids, GLboolean enabled)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glDebugMessageControlARB', source, type, severity, count, uint32(ids), enabled );

return
