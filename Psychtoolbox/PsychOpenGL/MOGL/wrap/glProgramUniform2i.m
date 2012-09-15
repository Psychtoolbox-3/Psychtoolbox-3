function glProgramUniform2i( program, location, v0, v1 )

% glProgramUniform2i  Interface to OpenGL function glProgramUniform2i
%
% usage:  glProgramUniform2i( program, location, v0, v1 )
%
% C function:  void glProgramUniform2i(GLuint program, GLint location, GLint v0, GLint v1)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2i', program, location, v0, v1 );

return
