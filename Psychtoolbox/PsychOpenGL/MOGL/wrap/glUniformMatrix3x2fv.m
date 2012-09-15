function glUniformMatrix3x2fv( location, count, transpose, value )

% glUniformMatrix3x2fv  Interface to OpenGL function glUniformMatrix3x2fv
%
% usage:  glUniformMatrix3x2fv( location, count, transpose, value )
%
% C function:  void glUniformMatrix3x2fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniformMatrix3x2fv', location, count, transpose, single(value) );

return
