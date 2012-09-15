function glVertexAttribP3uiv( index, type, normalized, value )

% glVertexAttribP3uiv  Interface to OpenGL function glVertexAttribP3uiv
%
% usage:  glVertexAttribP3uiv( index, type, normalized, value )
%
% C function:  void glVertexAttribP3uiv(GLuint index, GLenum type, GLboolean normalized, const GLuint* value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP3uiv', index, type, normalized, uint32(value) );

return
