function glVertexP4uiv( type, value )

% glVertexP4uiv  Interface to OpenGL function glVertexP4uiv
%
% usage:  glVertexP4uiv( type, value )
%
% C function:  void glVertexP4uiv(GLenum type, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP4uiv', type, uint32(value) );

return
