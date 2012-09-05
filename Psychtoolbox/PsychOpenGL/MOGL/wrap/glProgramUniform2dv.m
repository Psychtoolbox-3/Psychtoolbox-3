function glProgramUniform2dv( program, location, count, value )

% glProgramUniform2dv  Interface to OpenGL function glProgramUniform2dv
%
% usage:  glProgramUniform2dv( program, location, count, value )
%
% C function:  void glProgramUniform2dv(GLuint program, GLint location, GLsizei count, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniform2dv', program, location, count, double(value) );

return
