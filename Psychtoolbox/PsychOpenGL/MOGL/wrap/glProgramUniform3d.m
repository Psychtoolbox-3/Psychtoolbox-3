function glProgramUniform3d( program, location, v0, v1, v2 )

% glProgramUniform3d  Interface to OpenGL function glProgramUniform3d
%
% usage:  glProgramUniform3d( program, location, v0, v1, v2 )
%
% C function:  void glProgramUniform3d(GLuint program, GLint location, GLdouble v0, GLdouble v1, GLdouble v2)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform3d', program, location, v0, v1, v2 );

return
