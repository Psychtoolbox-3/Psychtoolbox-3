function glProgramUniform3uiv( program, location, count, value )

% glProgramUniform3uiv  Interface to OpenGL function glProgramUniform3uiv
%
% usage:  glProgramUniform3uiv( program, location, count, value )
%
% C function:  void glProgramUniform3uiv(GLuint program, GLint location, GLsizei count, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3uiv', program, location, count, uint32(value) );

return
