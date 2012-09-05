function glProgramUniform1ui( program, location, v0 )

% glProgramUniform1ui  Interface to OpenGL function glProgramUniform1ui
%
% usage:  glProgramUniform1ui( program, location, v0 )
%
% C function:  void glProgramUniform1ui(GLuint program, GLint location, GLuint v0)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1ui', program, location, v0 );

return
