function glProgramUniform3f( program, location, v0, v1, v2 )

% glProgramUniform3f  Interface to OpenGL function glProgramUniform3f
%
% usage:  glProgramUniform3f( program, location, v0, v1, v2 )
%
% C function:  void glProgramUniform3f(GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3f', program, location, v0, v1, v2 );

return
