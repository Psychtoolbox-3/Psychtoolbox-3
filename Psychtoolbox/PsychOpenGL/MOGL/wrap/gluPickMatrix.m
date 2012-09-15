function gluPickMatrix( x, y, delX, delY, viewport )

% gluPickMatrix  Interface to OpenGL function gluPickMatrix
%
% usage:  gluPickMatrix( x, y, delX, delY, viewport )
%
% C function:  void gluPickMatrix(GLdouble x, GLdouble y, GLdouble delX, GLdouble delY, GLint* viewport)

% 26-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'gluPickMatrix', x, y, delX, delY, viewport );

return
