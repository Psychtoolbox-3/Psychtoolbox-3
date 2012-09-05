function glProgramUniform4d( program, location, v0, v1, v2, v3 )

% glProgramUniform4d  Interface to OpenGL function glProgramUniform4d
%
% usage:  glProgramUniform4d( program, location, v0, v1, v2, v3 )
%
% C function:  void glProgramUniform4d(GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2, GLdouble v3)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4d', program, location, v0, v1, v2, v3 );

return
