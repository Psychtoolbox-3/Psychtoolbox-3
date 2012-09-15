function glUniformMatrix2x4fv( location, count, transpose, value )

% glUniformMatrix2x4fv  Interface to OpenGL function glUniformMatrix2x4fv
%
% usage:  glUniformMatrix2x4fv( location, count, transpose, value )
%
% C function:  void glUniformMatrix2x4fv(GLint location, GLsizei count, GLboolean transpose, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniformMatrix2x4fv', location, count, transpose, single(value) );

return
