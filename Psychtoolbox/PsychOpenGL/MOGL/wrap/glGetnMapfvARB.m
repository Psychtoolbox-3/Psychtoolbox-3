function glGetnMapfvARB( target, query, bufSize, v )


% glGetnMapfvARB  Interface to OpenGL function glGetnMapfvARB
%
% usage:  glGetnMapfvARB( target, query, bufSize, v )
%
% C function:  void glGetnMapfvARB(GLenum target, GLenum query, GLsizei bufSize, GLfloat* v)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glGetnMapfvARB', target, query, bufSize, v );

return
