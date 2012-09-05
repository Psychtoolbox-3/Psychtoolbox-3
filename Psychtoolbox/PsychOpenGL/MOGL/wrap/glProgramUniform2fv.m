function glProgramUniform2fv( program, location, count, value )

% glProgramUniform2fv  Interface to OpenGL function glProgramUniform2fv
%
% usage:  glProgramUniform2fv( program, location, count, value )
%
% C function:  void glProgramUniform2fv(GLuint program, GLint location, GLsizei count, const GLfloat* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2fv', program, location, count, single(value) );

return
