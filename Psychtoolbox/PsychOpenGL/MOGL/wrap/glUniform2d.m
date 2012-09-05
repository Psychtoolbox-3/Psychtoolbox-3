function glUniform2d( location, x, y )

% glUniform2d  Interface to OpenGL function glUniform2d
%
% usage:  glUniform2d( location, x, y )
%
% C function:  void glUniform2d(GLint location, GLdouble x, GLdouble y)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUniform2d', location, x, y );

return
