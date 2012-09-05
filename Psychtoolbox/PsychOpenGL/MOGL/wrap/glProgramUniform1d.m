function glProgramUniform1d( program, location, v0 )

% glProgramUniform1d  Interface to OpenGL function glProgramUniform1d
%
% usage:  glProgramUniform1d( program, location, v0 )
%
% C function:  void glProgramUniform1d(GLuint program, GLint location, GLdouble v0)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1d', program, location, v0 );

return
