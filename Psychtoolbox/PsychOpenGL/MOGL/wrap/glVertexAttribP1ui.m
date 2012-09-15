function glVertexAttribP1ui( index, type, normalized, value )

% glVertexAttribP1ui  Interface to OpenGL function glVertexAttribP1ui
%
% usage:  glVertexAttribP1ui( index, type, normalized, value )
%
% C function:  void glVertexAttribP1ui(GLuint index, GLenum type, GLboolean normalized, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribP1ui', index, type, normalized, value );

return
