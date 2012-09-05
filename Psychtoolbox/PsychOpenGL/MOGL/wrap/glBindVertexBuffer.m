function glBindVertexBuffer( bindingindex, buffer, ptr, stride )

% glBindVertexBuffer  Interface to OpenGL function glBindVertexBuffer
%
% usage:  glBindVertexBuffer( bindingindex, buffer, ptr, stride )
%
% C function:  void glBindVertexBuffer(GLuint bindingindex, GLuint buffer, GLint ptr, GLsizei stride)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glBindVertexBuffer', bindingindex, buffer, ptr, stride );

return
