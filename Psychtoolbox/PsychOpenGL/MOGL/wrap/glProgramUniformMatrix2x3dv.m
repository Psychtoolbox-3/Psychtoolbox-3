function glProgramUniformMatrix2x3dv( program, location, count, transpose, value )

% glProgramUniformMatrix2x3dv  Interface to OpenGL function glProgramUniformMatrix2x3dv
%
% usage:  glProgramUniformMatrix2x3dv( program, location, count, transpose, value )
%
% C function:  void glProgramUniformMatrix2x3dv(GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniformMatrix2x3dv', program, location, count, transpose, double(value) );

return
