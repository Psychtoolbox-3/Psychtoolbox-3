function glVertexAttribP3ui( index, type, normalized, value )

% glVertexAttribP3ui  Interface to OpenGL function glVertexAttribP3ui
%
% usage:  glVertexAttribP3ui( index, type, normalized, value )
%
% C function:  void glVertexAttribP3ui(GLuint index, GLenum type, GLboolean normalized, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP3ui', index, type, normalized, value );

return
