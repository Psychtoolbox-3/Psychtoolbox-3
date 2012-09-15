function glProgramUniform1fv( program, location, count, value )

% glProgramUniform1fv  Interface to OpenGL function glProgramUniform1fv
%
% usage:  glProgramUniform1fv( program, location, count, value )
%
% C function:  void glProgramUniform1fv(GLuint program, GLint location, GLsizei count, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1fv', program, location, count, single(value) );

return
