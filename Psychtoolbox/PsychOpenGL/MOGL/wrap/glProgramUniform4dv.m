function glProgramUniform4dv( program, location, count, value )

% glProgramUniform4dv  Interface to OpenGL function glProgramUniform4dv
%
% usage:  glProgramUniform4dv( program, location, count, value )
%
% C function:  void glProgramUniform4dv(GLuint program, GLint location, GLsizei count, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform4dv', program, location, count, double(value) );

return
