function glOrthof( left, right, bottom, top, zNear, zFar )

% glOrthof  Interface to OpenGL-ES function glOrthof
%
% usage:  glOrthof( left, right, bottom, top, zNear, zFar )
%
% C function:  void glOrthof(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat zNear, GLfloat zFar)

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

m = zeros(4,4);
m(1,1) = 2 / (right - left);
m(2,2) = 2 / (top - bottom);
m(3,3) = -2 / (zFar - zNear);
m(4,4) = 1;
m(1,4) = -(right + left) / (right - left);
m(2,4) = -(top + bottom) / (top - bottom);
m(3,4) = -(zFar + zNear) / (zFar - zNear);

% TODO FIXME: Need m transposed m' ?
glMultMatrixf(m);

return
