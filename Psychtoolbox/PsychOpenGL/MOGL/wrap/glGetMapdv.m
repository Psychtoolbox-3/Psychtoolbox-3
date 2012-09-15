function glGetMapdv( target, query, v )

% glGetMapdv  Interface to OpenGL function glGetMapdv
%
% usage:  glGetMapdv( target, query, v )
%
% C function:  void glGetMapdv(GLenum target, GLenum query, GLdouble* v)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glGetMapdv', target, query, v );

return
