function glVertexAttribP4ui( index, type, normalized, value )

% glVertexAttribP4ui  Interface to OpenGL function glVertexAttribP4ui
%
% usage:  glVertexAttribP4ui( index, type, normalized, value )
%
% C function:  void glVertexAttribP4ui(GLuint index, GLenum type, GLboolean normalized, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP4ui', index, type, normalized, value );

return
