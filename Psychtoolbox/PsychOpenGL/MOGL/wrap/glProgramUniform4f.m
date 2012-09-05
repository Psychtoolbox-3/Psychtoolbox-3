function glProgramUniform4f( program, location, v0, v1, v2, v3 )

% glProgramUniform4f  Interface to OpenGL function glProgramUniform4f
%
% usage:  glProgramUniform4f( program, location, v0, v1, v2, v3 )
%
% C function:  void glProgramUniform4f(GLuint program, GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4f', program, location, v0, v1, v2, v3 );

return
