function gluOrtho2D( left, right, bottom, top )

% gluOrtho2D  Interface to OpenGL function gluOrtho2D
%
% usage:  gluOrtho2D( left, right, bottom, top )
%
% C function:  void gluOrtho2D(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top)

% 26-Mar-2011 -- created (generated automatically from header files)
%  4-Apr-2013 Made compatible with OpenGL-ES.

% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'gluOrtho2D', left, right, bottom, top );
    return;
end

glOrthof(left, right, bottom, top, -1, +1);

return
