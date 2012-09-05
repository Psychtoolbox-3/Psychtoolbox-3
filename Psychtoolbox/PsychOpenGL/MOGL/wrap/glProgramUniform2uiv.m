function glProgramUniform2uiv( program, location, count, value )

% glProgramUniform2uiv  Interface to OpenGL function glProgramUniform2uiv
%
% usage:  glProgramUniform2uiv( program, location, count, value )
%
% C function:  void glProgramUniform2uiv(GLuint program, GLint location, GLsizei count, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2uiv', program, location, count, uint32(value) );

return
