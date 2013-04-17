function glFrustum( left, right, bottom, top, zNear, zFar )

% glFrustum  Interface to OpenGL function glFrustum
%
% usage:  glFrustum( left, right, bottom, top, zNear, zFar )
%
% C function:  void glFrustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar)

% 25-Mar-2011 -- created (generated automatically from header files)
%  4-Apr-2013 Made compatible with OpenGL-ES.

% ---protected---

if nargin~=6,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glFrustum', left, right, bottom, top, zNear, zFar );
    return;
end

% Manual definition:
m = zeros(4,4);
m(1,1) = 2 * zNear / (right - left);
m(2,2) = 2 * zNear / (top - bottom);
m(1,3) = (right + left) / (right - left);
m(2,3) = (top + bottom) / (top - bottom);
m(3,3) = - (zFar + zNear) / (zFar - zNear);
m(3,4) = - 2 * zFar * zNear / (zFar - zNear);
m(4,3) = -1;

% TODO FIXME: Check if this must be m or m' aka m transposed?
glMultMatrixf(m);

return
