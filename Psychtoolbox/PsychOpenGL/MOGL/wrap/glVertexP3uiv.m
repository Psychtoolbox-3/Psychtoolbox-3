function glVertexP3uiv( type, value )

% glVertexP3uiv  Interface to OpenGL function glVertexP3uiv
%
% usage:  glVertexP3uiv( type, value )
%
% C function:  void glVertexP3uiv(GLenum type, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP3uiv', type, uint32(value) );

return
