function glProgramUniform3ui( program, location, v0, v1, v2 )

% glProgramUniform3ui  Interface to OpenGL function glProgramUniform3ui
%
% usage:  glProgramUniform3ui( program, location, v0, v1, v2 )
%
% C function:  void glProgramUniform3ui(GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3ui', program, location, v0, v1, v2 );

return
