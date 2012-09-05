function glProgramUniform2iv( program, location, count, value )

% glProgramUniform2iv  Interface to OpenGL function glProgramUniform2iv
%
% usage:  glProgramUniform2iv( program, location, count, value )
%
% C function:  void glProgramUniform2iv(GLuint program, GLint location, GLsizei count, const GLint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2iv', program, location, count, int32(value) );

return
