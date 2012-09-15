function glProgramUniform4iv( program, location, count, value )

% glProgramUniform4iv  Interface to OpenGL function glProgramUniform4iv
%
% usage:  glProgramUniform4iv( program, location, count, value )
%
% C function:  void glProgramUniform4iv(GLuint program, GLint location, GLsizei count, const GLint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4iv', program, location, count, int32(value) );

return
