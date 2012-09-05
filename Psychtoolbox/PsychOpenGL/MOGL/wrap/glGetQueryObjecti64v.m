function params = glGetQueryObjecti64v( id, pname )

% glGetQueryObjecti64v  Interface to OpenGL function glGetQueryObjecti64v
%
% usage:  params = glGetQueryObjecti64v( id, pname )
%
% C function:  void glGetQueryObjecti64v(GLuint id, GLenum pname, GLint64* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int64(0);

moglcore( 'glGetQueryObjecti64v', id, pname, params );

return
