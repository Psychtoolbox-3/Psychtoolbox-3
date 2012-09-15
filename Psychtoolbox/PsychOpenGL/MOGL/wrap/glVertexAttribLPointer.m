function glVertexAttribLPointer( index, size, type, stride, pointer )

% glVertexAttribLPointer  Interface to OpenGL function glVertexAttribLPointer
%
% usage:  glVertexAttribLPointer( index, size, type, stride, pointer )
%
% C function:  void glVertexAttribLPointer(GLuint index, GLint size, GLenum type, GLsizei stride, const GLvoid* pointer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribLPointer', index, size, type, stride, pointer );

return
