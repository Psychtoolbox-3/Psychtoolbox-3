function glUniform3d( location, x, y, z )

% glUniform3d  Interface to OpenGL function glUniform3d
%
% usage:  glUniform3d( location, x, y, z )
%
% C function:  void glUniform3d(GLint location, GLdouble x, GLdouble y, GLdouble z)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniform3d', location, x, y, z );

return
