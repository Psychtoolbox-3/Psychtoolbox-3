function glMultiTexCoordP2ui( texture, type, coords )

% glMultiTexCoordP2ui  Interface to OpenGL function glMultiTexCoordP2ui
%
% usage:  glMultiTexCoordP2ui( texture, type, coords )
%
% C function:  void glMultiTexCoordP2ui(GLenum texture, GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP2ui', texture, type, coords );

return
