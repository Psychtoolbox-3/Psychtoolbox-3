function glProgramUniform4i( program, location, v0, v1, v2, v3 )

% glProgramUniform4i  Interface to OpenGL function glProgramUniform4i
%
% usage:  glProgramUniform4i( program, location, v0, v1, v2, v3 )
%
% C function:  void glProgramUniform4i(GLuint program, GLint location, GLint v0, GLint v1, GLint v2, GLint v3)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4i', program, location, v0, v1, v2, v3 );

return
