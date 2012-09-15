function glProgramUniform1dv( program, location, count, value )

% glProgramUniform1dv  Interface to OpenGL function glProgramUniform1dv
%
% usage:  glProgramUniform1dv( program, location, count, value )
%
% C function:  void glProgramUniform1dv(GLuint program, GLint location, GLsizei count, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform1dv', program, location, count, double(value) );

return
