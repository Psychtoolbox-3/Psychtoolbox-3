function glProgramUniform3fv( program, location, count, value )

% glProgramUniform3fv  Interface to OpenGL function glProgramUniform3fv
%
% usage:  glProgramUniform3fv( program, location, count, value )
%
% C function:  void glProgramUniform3fv(GLuint program, GLint location, GLsizei count, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3fv', program, location, count, single(value) );

return
