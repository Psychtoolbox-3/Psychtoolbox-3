function glProgramUniform3iv( program, location, count, value )

% glProgramUniform3iv  Interface to OpenGL function glProgramUniform3iv
%
% usage:  glProgramUniform3iv( program, location, count, value )
%
% C function:  void glProgramUniform3iv(GLuint program, GLint location, GLsizei count, const GLint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3iv', program, location, count, int32(value) );

return
