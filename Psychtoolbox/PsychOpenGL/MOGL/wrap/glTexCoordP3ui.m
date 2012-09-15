function glTexCoordP3ui( type, coords )

% glTexCoordP3ui  Interface to OpenGL function glTexCoordP3ui
%
% usage:  glTexCoordP3ui( type, coords )
%
% C function:  void glTexCoordP3ui(GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP3ui', type, coords );

return
