function glGetMapiv( target, query, v )

% glGetMapiv  Interface to OpenGL function glGetMapiv
%
% usage:  glGetMapiv( target, query, v )
%
% C function:  void glGetMapiv(GLenum target, GLenum query, GLint* v)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glGetMapiv', target, query, v );

return
