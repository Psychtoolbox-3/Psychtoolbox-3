function glUniform1d( location, x )

% glUniform1d  Interface to OpenGL function glUniform1d
%
% usage:  glUniform1d( location, x )
%
% C function:  void glUniform1d(GLint location, GLdouble x)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glUniform1d', location, x );

return
