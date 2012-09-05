function glViewportArrayv( first, count, v )

% glViewportArrayv  Interface to OpenGL function glViewportArrayv
%
% usage:  glViewportArrayv( first, count, v )
%
% C function:  void glViewportArrayv(GLuint first, GLsizei count, const GLfloat* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glViewportArrayv', first, count, single(v) );

return
