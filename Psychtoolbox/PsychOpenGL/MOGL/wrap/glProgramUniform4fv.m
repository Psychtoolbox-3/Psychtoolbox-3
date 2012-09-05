function glProgramUniform4fv( program, location, count, value )

% glProgramUniform4fv  Interface to OpenGL function glProgramUniform4fv
%
% usage:  glProgramUniform4fv( program, location, count, value )
%
% C function:  void glProgramUniform4fv(GLuint program, GLint location, GLsizei count, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4fv', program, location, count, single(value) );

return
