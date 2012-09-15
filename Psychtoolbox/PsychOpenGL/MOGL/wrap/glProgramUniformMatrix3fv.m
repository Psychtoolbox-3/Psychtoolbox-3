function glProgramUniformMatrix3fv( program, location, count, transpose, value )

% glProgramUniformMatrix3fv  Interface to OpenGL function glProgramUniformMatrix3fv
%
% usage:  glProgramUniformMatrix3fv( program, location, count, transpose, value )
%
% C function:  void glProgramUniformMatrix3fv(GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniformMatrix3fv', program, location, count, transpose, single(value) );

return
