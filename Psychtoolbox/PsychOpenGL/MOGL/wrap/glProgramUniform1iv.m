function glProgramUniform1iv( program, location, count, value )

% glProgramUniform1iv  Interface to OpenGL function glProgramUniform1iv
%
% usage:  glProgramUniform1iv( program, location, count, value )
%
% C function:  void glProgramUniform1iv(GLuint program, GLint location, GLsizei count, const GLint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1iv', program, location, count, int32(value) );

return
