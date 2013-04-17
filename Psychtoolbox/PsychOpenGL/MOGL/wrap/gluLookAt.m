function gluLookAt( eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ )

% gluLookAt  Interface to OpenGL function gluLookAt
%
% usage:  gluLookAt( eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ )
%
% C function:  void gluLookAt(GLdouble eyeX, GLdouble eyeY, GLdouble eyeZ, GLdouble centerX, GLdouble centerY, GLdouble centerZ, GLdouble upX, GLdouble upY, GLdouble upZ)

% 25-Mar-2011 -- created (generated automatically from header files)
%  4-Apr-2013 -- Manually added Matlab implementation for OpenGL-ES1.x

% ---protected---

if nargin~=9,
    error('invalid number of arguments');
end

% OpenGL-ES1 active?
if ~IsGLES1
    % No. Dispatch to GLU:
    moglcore( 'gluLookAt', eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ );
    return;
end

% Yes. Need to implement gluLookAt() manually here.
forward(1) = centerX - eyeX;
forward(2) = centerY - eyeY;
forward(3) = centerZ - eyeZ;
forward = normalize(forward);

up(1) = upX;
up(2) = upY;
up(3) = upZ;
up = normalize(up);

% Side = forward x up
side = cross(forward, up);
side = normalize(side);

% Recompute up as: up = side x forward
up = cross(side, forward);

m = diag([1 1 1 1]);
m(1,1) = side(1);
m(2,1) = side(2);
m(3,1) = side(3);

m(1,2) = up(1);
m(2,2) = up(2);
m(3,2) = up(3);

m(1,3) = -forward(1);
m(2,3) = -forward(2);
m(3,3) = -forward(3);

glMultMatrixf(m');
glTranslatef(-eyeX, -eyeY, -eyeZ);

return;

% Normalize vector to unit-length:
function v = normalize(v)
    l = norm(v);
    if l == 0
        return;
    end

    v = v / l;
return;
