function glVertexP2uiv( type, value )

% glVertexP2uiv  Interface to OpenGL function glVertexP2uiv
%
% usage:  glVertexP2uiv( type, value )
%
% C function:  void glVertexP2uiv(GLenum type, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP2uiv', type, uint32(value) );

return
