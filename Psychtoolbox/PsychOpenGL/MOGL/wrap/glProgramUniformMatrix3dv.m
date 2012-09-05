function glProgramUniformMatrix3dv( program, location, count, transpose, value )

% glProgramUniformMatrix3dv  Interface to OpenGL function glProgramUniformMatrix3dv
%
% usage:  glProgramUniformMatrix3dv( program, location, count, transpose, value )
%
% C function:  void glProgramUniformMatrix3dv(GLuint program, GLint location, GLsizei count, GLboolean transpose, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glProgramUniformMatrix3dv', program, location, count, transpose, double(value) );

return
