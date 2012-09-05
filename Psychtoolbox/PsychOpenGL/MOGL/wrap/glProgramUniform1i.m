function glProgramUniform1i( program, location, v0 )

% glProgramUniform1i  Interface to OpenGL function glProgramUniform1i
%
% usage:  glProgramUniform1i( program, location, v0 )
%
% C function:  void glProgramUniform1i(GLuint program, GLint location, GLint v0)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1i', program, location, v0 );

return
