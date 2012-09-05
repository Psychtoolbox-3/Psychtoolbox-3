function glProgramUniform2ui( program, location, v0, v1 )

% glProgramUniform2ui  Interface to OpenGL function glProgramUniform2ui
%
% usage:  glProgramUniform2ui( program, location, v0, v1 )
%
% C function:  void glProgramUniform2ui(GLuint program, GLint location, GLuint v0, GLuint v1)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2ui', program, location, v0, v1 );

return
