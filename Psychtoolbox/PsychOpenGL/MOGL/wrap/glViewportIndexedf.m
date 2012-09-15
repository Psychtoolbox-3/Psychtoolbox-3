function glViewportIndexedf( index, x, y, w, h )

% glViewportIndexedf  Interface to OpenGL function glViewportIndexedf
%
% usage:  glViewportIndexedf( index, x, y, w, h )
%
% C function:  void glViewportIndexedf(GLuint index, GLfloat x, GLfloat y, GLfloat w, GLfloat h)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glViewportIndexedf', index, x, y, w, h );

return
