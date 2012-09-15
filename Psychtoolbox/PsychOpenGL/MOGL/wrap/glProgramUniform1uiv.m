function glProgramUniform1uiv( program, location, count, value )

% glProgramUniform1uiv  Interface to OpenGL function glProgramUniform1uiv
%
% usage:  glProgramUniform1uiv( program, location, count, value )
%
% C function:  void glProgramUniform1uiv(GLuint program, GLint location, GLsizei count, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1uiv', program, location, count, uint32(value) );

return
