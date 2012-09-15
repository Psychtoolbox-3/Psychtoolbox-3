function val = glGetMultisamplefv( pname, index )

% glGetMultisamplefv  Interface to OpenGL function glGetMultisamplefv
%
% usage:  val = glGetMultisamplefv( pname, index )
%
% C function:  void glGetMultisamplefv(GLenum pname, GLuint index, GLfloat* val)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

val = single([0,0]);

moglcore( 'glGetMultisamplefv', pname, index, val );

return
