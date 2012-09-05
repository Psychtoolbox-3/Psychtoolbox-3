function glVertexAttribP4uiv( index, type, normalized, value )

% glVertexAttribP4uiv  Interface to OpenGL function glVertexAttribP4uiv
%
% usage:  glVertexAttribP4uiv( index, type, normalized, value )
%
% C function:  void glVertexAttribP4uiv(GLuint index, GLenum type, GLboolean normalized, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP4uiv', index, type, normalized, uint32(value) );

return
