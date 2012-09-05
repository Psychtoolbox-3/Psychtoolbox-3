function glMultiTexCoordP2uiv( texture, type, coords )

% glMultiTexCoordP2uiv  Interface to OpenGL function glMultiTexCoordP2uiv
%
% usage:  glMultiTexCoordP2uiv( texture, type, coords )
%
% C function:  void glMultiTexCoordP2uiv(GLenum texture, GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP2uiv', texture, type, uint32(coords) );

return
