function glTexCoordP4ui( type, coords )

% glTexCoordP4ui  Interface to OpenGL function glTexCoordP4ui
%
% usage:  glTexCoordP4ui( type, coords )
%
% C function:  void glTexCoordP4ui(GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP4ui', type, coords );

return
