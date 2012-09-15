function glSecondaryColorP3uiv( type, color )

% glSecondaryColorP3uiv  Interface to OpenGL function glSecondaryColorP3uiv
%
% usage:  glSecondaryColorP3uiv( type, color )
%
% C function:  void glSecondaryColorP3uiv(GLenum type, const GLuint* color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glSecondaryColorP3uiv', type, uint32(color) );

return
