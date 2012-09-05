function glTexCoordP1uiv( type, coords )

% glTexCoordP1uiv  Interface to OpenGL function glTexCoordP1uiv
%
% usage:  glTexCoordP1uiv( type, coords )
%
% C function:  void glTexCoordP1uiv(GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP1uiv', type, uint32(coords) );

return
