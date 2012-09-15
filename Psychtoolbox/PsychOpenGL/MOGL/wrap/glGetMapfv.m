function glGetMapfv( target, query, v )

% glGetMapfv  Interface to OpenGL function glGetMapfv
%
% usage: glGetMapfv( target, query, v )
%
% C function:  void glGetMapfv(GLenum target, GLenum query, GLfloat* v)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glGetMapfv', target, query, v );

return
