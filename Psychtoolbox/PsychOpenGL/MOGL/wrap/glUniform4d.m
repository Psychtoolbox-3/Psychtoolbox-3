function glUniform4d( location, x, y, z, w )

% glUniform4d  Interface to OpenGL function glUniform4d
%
% usage:  glUniform4d( location, x, y, z, w )
%
% C function:  void glUniform4d(GLint location, GLdouble x, GLdouble y, GLdouble z, GLdouble w)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glUniform4d', location, x, y, z, w );

return
