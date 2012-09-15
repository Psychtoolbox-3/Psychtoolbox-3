function glUniform3dv( location, count, value )

% glUniform3dv  Interface to OpenGL function glUniform3dv
%
% usage:  glUniform3dv( location, count, value )
%
% C function:  void glUniform3dv(GLint location, GLsizei count, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUniform3dv', location, count, double(value) );

return
