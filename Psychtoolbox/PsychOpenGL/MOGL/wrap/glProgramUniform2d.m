function glProgramUniform2d( program, location, v0, v1 )

% glProgramUniform2d  Interface to OpenGL function glProgramUniform2d
%
% usage:  glProgramUniform2d( program, location, v0, v1 )
%
% C function:  void glProgramUniform2d(GLuint program, GLint location, GLdouble v0, GLdouble v1)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2d', program, location, v0, v1 );

return
