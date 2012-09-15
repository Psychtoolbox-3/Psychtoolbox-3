function params = glGetQueryObjectui64v( id, pname )

% glGetQueryObjectui64v  Interface to OpenGL function glGetQueryObjectui64v
%
% usage:  params = glGetQueryObjectui64v( id, pname )
%
% C function:  void glGetQueryObjectui64v(GLuint id, GLenum pname, GLuint64* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = uint64(0);

moglcore( 'glGetQueryObjectui64v', id, pname, params );

return
