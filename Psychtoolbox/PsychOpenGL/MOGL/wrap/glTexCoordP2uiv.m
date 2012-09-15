function glTexCoordP2uiv( type, coords )

% glTexCoordP2uiv  Interface to OpenGL function glTexCoordP2uiv
%
% usage:  glTexCoordP2uiv( type, coords )
%
% C function:  void glTexCoordP2uiv(GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP2uiv', type, uint32(coords) );

return
