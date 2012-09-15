function glVertexAttribP2ui( index, type, normalized, value )

% glVertexAttribP2ui  Interface to OpenGL function glVertexAttribP2ui
%
% usage:  glVertexAttribP2ui( index, type, normalized, value )
%
% C function:  void glVertexAttribP2ui(GLuint index, GLenum type, GLboolean normalized, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP2ui', index, type, normalized, value );

return
