function glMultiTexCoordP4uiv( texture, type, coords )

% glMultiTexCoordP4uiv  Interface to OpenGL function glMultiTexCoordP4uiv
%
% usage:  glMultiTexCoordP4uiv( texture, type, coords )
%
% C function:  void glMultiTexCoordP4uiv(GLenum texture, GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP4uiv', texture, type, uint32(coords) );

return
