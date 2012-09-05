function glGetnMapdvARB( target, query, bufSize, v )

% glGetnMapdvARB  Interface to OpenGL function glGetnMapdvARB
%
% usage:  glGetnMapdvARB( target, query, bufSize, v )
%
% C function:  void glGetnMapdvARB(GLenum target, GLenum query, GLsizei bufSize, GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glGetnMapdvARB', target, query, bufSize, v );

return
