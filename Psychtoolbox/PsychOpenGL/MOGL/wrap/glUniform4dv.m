function glUniform4dv( location, count, value )

% glUniform4dv  Interface to OpenGL function glUniform4dv
%
% usage:  glUniform4dv( location, count, value )
%
% C function:  void glUniform4dv(GLint location, GLsizei count, const GLdouble* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUniform4dv', location, count, double(value) );

return
