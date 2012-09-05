function glUniformMatrix2x3dv( location, count, transpose, value )

% glUniformMatrix2x3dv  Interface to OpenGL function glUniformMatrix2x3dv
%
% usage:  glUniformMatrix2x3dv( location, count, transpose, value )
%
% C function:  void glUniformMatrix2x3dv(GLint location, GLsizei count, GLboolean transpose, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glUniformMatrix2x3dv', location, count, transpose, double(value) );

return
