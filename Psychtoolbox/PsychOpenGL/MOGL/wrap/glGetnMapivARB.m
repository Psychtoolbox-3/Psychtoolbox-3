function glGetnMapivARB( target, query, bufSize, v )
% glGetnMapivARB  Interface to OpenGL function glGetnMapivARB
%
% usage:  glGetnMapivARB( target, query, bufSize, v )
%
% C function:  void glGetnMapivARB(GLenum target, GLenum query, GLsizei bufSize, GLint* v)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glGetnMapivARB', target, query, bufSize, v );

return
