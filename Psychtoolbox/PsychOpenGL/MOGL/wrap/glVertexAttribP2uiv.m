function glVertexAttribP2uiv( index, type, normalized, value )

% glVertexAttribP2uiv  Interface to OpenGL function glVertexAttribP2uiv
%
% usage:  glVertexAttribP2uiv( index, type, normalized, value )
%
% C function:  void glVertexAttribP2uiv(GLuint index, GLenum type, GLboolean normalized, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP2uiv', index, type, normalized, uint32(value) );

return
