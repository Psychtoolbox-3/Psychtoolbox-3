function glMultiTexCoordP3uiv( texture, type, coords )

% glMultiTexCoordP3uiv  Interface to OpenGL function glMultiTexCoordP3uiv
%
% usage:  glMultiTexCoordP3uiv( texture, type, coords )
%
% C function:  void glMultiTexCoordP3uiv(GLenum texture, GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glMultiTexCoordP3uiv', texture, type, uint32(coords) );

return
