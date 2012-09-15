function glProgramUniformMatrix3x4fv( program, location, count, transpose, value )

% glProgramUniformMatrix3x4fv  Interface to OpenGL function glProgramUniformMatrix3x4fv
%
% usage:  glProgramUniformMatrix3x4fv( program, location, count, transpose, value )
%
% C function:  void glProgramUniformMatrix3x4fv(GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniformMatrix3x4fv', program, location, count, transpose, single(value) );

return
