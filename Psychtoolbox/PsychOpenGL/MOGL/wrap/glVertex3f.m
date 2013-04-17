function glVertex3f( x, y, z )

% glVertex3f  Interface to OpenGL function glVertex3f
%
% usage:  glVertex3f( x, y, z )
%
% C function:  void glVertex3f(GLfloat x, GLfloat y, GLfloat z)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glVertex3f', x, y, z );
else
    moglcore( 'ftglVertex3f', x, y, z );
end

return
