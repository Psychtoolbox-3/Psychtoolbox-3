function glUniformMatrix3dv( location, count, transpose, value )

% glUniformMatrix3dv  Interface to OpenGL function glUniformMatrix3dv
%
% usage:  glUniformMatrix3dv( location, count, transpose, value )
%
% C function:  void glUniformMatrix3dv(GLint location, GLsizei count, GLboolean transpose, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniformMatrix3dv', location, count, transpose, double(value) );

return
