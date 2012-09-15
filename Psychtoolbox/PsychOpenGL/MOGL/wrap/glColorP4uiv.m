function glColorP4uiv( type, color )

% glColorP4uiv  Interface to OpenGL function glColorP4uiv
%
% usage:  glColorP4uiv( type, color )
%
% C function:  void glColorP4uiv(GLenum type, const GLuint* color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glColorP4uiv', type, uint32(color) );

return
