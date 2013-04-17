function glOrtho( left, right, bottom, top, zNear, zFar )

% glOrtho  Interface to OpenGL function glOrtho
%
% usage:  glOrtho( left, right, bottom, top, zNear, zFar )
%
% C function:  void glOrtho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar)

% 25-Mar-2011 -- created (generated automatically from header files)
%  4-Apr-2013 Made compatible with OpenGL-ES.

% ---protected---

if nargin~=6,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glOrtho', left, right, bottom, top, zNear, zFar );
    return;
end

glOrthof(left, right, bottom, top, zNear, zFar);

return
