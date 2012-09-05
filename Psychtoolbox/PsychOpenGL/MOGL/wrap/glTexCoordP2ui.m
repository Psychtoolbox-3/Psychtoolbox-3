function glTexCoordP2ui( type, coords )

% glTexCoordP2ui  Interface to OpenGL function glTexCoordP2ui
%
% usage:  glTexCoordP2ui( type, coords )
%
% C function:  void glTexCoordP2ui(GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP2ui', type, coords );

return
