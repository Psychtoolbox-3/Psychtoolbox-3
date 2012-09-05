function params = glGetQueryIndexediv( target, index, pname )

% glGetQueryIndexediv  Interface to OpenGL function glGetQueryIndexediv
%
% usage:  params = glGetQueryIndexediv( target, index, pname )
%
% C function:  void glGetQueryIndexediv(GLenum target, GLuint index, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetQueryIndexediv', target, index, pname, params );

return
