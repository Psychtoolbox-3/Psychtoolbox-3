function glProgramUniform1f( program, location, v0 )

% glProgramUniform1f  Interface to OpenGL function glProgramUniform1f
%
% usage:  glProgramUniform1f( program, location, v0 )
%
% C function:  void glProgramUniform1f(GLuint program, GLint location, GLfloat v0)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1f', program, location, v0 );

return
