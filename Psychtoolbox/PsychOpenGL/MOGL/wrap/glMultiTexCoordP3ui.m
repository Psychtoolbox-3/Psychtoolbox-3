function glMultiTexCoordP3ui( texture, type, coords )

% glMultiTexCoordP3ui  Interface to OpenGL function glMultiTexCoordP3ui
%
% usage:  glMultiTexCoordP3ui( texture, type, coords )
%
% C function:  void glMultiTexCoordP3ui(GLenum texture, GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP3ui', texture, type, coords );

return
