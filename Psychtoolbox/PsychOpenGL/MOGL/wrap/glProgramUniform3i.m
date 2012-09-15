function glProgramUniform3i( program, location, v0, v1, v2 )

% glProgramUniform3i  Interface to OpenGL function glProgramUniform3i
%
% usage:  glProgramUniform3i( program, location, v0, v1, v2 )
%
% C function:  void glProgramUniform3i(GLuint program, GLint location, GLint v0, GLint v1, GLint v2)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3i', program, location, v0, v1, v2 );

return
