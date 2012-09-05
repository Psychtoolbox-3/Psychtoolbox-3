function params = glGetQueryObjectiv( id, pname )

% glGetQueryObjectiv  Interface to OpenGL function glGetQueryObjectiv
%
% usage:  params = glGetQueryObjectiv( id, pname )
%
% C function:  void glGetQueryObjectiv(GLuint id, GLenum pname, GLint* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetQueryObjectiv', id, pname, params );

return
