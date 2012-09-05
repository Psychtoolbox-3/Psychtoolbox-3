function params = glGetQueryObjectuiv( id, pname )

% glGetQueryObjectuiv  Interface to OpenGL function glGetQueryObjectuiv
%
% usage:  params = glGetQueryObjectuiv( id, pname )
%
% C function:  void glGetQueryObjectuiv(GLuint id, GLenum pname, GLuint* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = uint32(0);

moglcore( 'glGetQueryObjectuiv', id, pname, params );

return
