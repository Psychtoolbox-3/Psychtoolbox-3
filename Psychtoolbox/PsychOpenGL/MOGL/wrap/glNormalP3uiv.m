function glNormalP3uiv( type, coords )

% glNormalP3uiv  Interface to OpenGL function glNormalP3uiv
%
% usage:  glNormalP3uiv( type, coords )
%
% C function:  void glNormalP3uiv(GLenum type, const GLuint* coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glNormalP3uiv', type, uint32(coords) );

return
