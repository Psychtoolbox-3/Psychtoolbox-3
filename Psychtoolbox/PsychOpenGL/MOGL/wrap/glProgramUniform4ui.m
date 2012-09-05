function glProgramUniform4ui( program, location, v0, v1, v2, v3 )

% glProgramUniform4ui  Interface to OpenGL function glProgramUniform4ui
%
% usage:  glProgramUniform4ui( program, location, v0, v1, v2, v3 )
%
% C function:  void glProgramUniform4ui(GLuint program, GLint location, GLuint v0, GLuint v1, GLuint v2, GLuint v3)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4ui', program, location, v0, v1, v2, v3 );

return
