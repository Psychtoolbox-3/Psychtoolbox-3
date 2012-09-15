function glMultiTexCoordP1uiv( texture, type, coords )

% glMultiTexCoordP1uiv  Interface to OpenGL function glMultiTexCoordP1uiv
%
% usage:  glMultiTexCoordP1uiv( texture, type, coords )
%
% C function:  void glMultiTexCoordP1uiv(GLenum texture, GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP1uiv', texture, type, uint32(coords) );

return
