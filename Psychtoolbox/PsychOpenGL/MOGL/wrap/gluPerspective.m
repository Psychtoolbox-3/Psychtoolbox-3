function gluPerspective( fovy, aspect, zNear, zFar )

% gluPerspective  Interface to OpenGL function gluPerspective
%
% usage:  gluPerspective( fovy, aspect, zNear, zFar )
%
% C function:  void gluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)

% 26-Mar-2011 -- created (generated automatically from header files)
%  4-Apr-2013 -- Manually added Matlab implementation for OpenGL-ES1.x

% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

% OpenGL-ES1 active?
if ~IsGLES1
    % No. Dispatch to GLU:
    moglcore( 'gluPerspective', fovy, aspect, zNear, zFar );
    return;
end

% Yes. Need to implement gluPerspective() manually here.
m = zeros(4,4);
f = cot(fovy * 2 * pi / 360 / 2);
m(1,1) = f / aspect;
m(2,2) = f;
m(3,3) = (zFar + zNear) / (zNear - zFar);
m(3,4) = 2 * zFar * zNear / (zNear - zFar);
m(4,3) = -1;
glMultMatrixf(m);

return;
