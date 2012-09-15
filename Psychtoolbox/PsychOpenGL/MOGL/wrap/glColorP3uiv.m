function glColorP3uiv( type, color )

% glColorP3uiv  Interface to OpenGL function glColorP3uiv
%
% usage:  glColorP3uiv( type, color )
%
% C function:  void glColorP3uiv(GLenum type, const GLuint* color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glColorP3uiv', type, uint32(color) );

return
