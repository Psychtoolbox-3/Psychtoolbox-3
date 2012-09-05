function glTexCoordP3uiv( type, coords )

% glTexCoordP3uiv  Interface to OpenGL function glTexCoordP3uiv
%
% usage:  glTexCoordP3uiv( type, coords )
%
% C function:  void glTexCoordP3uiv(GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP3uiv', type, uint32(coords) );

return
