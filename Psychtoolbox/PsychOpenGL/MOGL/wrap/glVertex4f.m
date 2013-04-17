function glVertex4f( x, y, z, w )

% glVertex4f  Interface to OpenGL function glVertex4f
%
% usage:  glVertex4f( x, y, z, w )
%
% C function:  void glVertex4f(GLfloat x, GLfloat y, GLfloat z, GLfloat w)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glVertex4f', x, y, z, w );
else
    moglcore( 'ftglVertex4f', x, y, z, w );
end

return
