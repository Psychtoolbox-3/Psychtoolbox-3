function glUniformMatrix3x2dv( location, count, transpose, value )

% glUniformMatrix3x2dv  Interface to OpenGL function glUniformMatrix3x2dv
%
% usage:  glUniformMatrix3x2dv( location, count, transpose, value )
%
% C function:  void glUniformMatrix3x2dv(GLint location, GLsizei count, GLboolean transpose, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniformMatrix3x2dv', location, count, transpose, double(value) );

return
