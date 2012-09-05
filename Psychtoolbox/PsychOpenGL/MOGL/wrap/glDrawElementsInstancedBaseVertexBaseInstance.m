function glDrawElementsInstancedBaseVertexBaseInstance( mode, count, type, indices, instancecount, basevertex, baseinstance )

% glDrawElementsInstancedBaseVertexBaseInstance  Interface to OpenGL function glDrawElementsInstancedBaseVertexBaseInstance
%
% usage:  glDrawElementsInstancedBaseVertexBaseInstance( mode, count, type, indices, instancecount, basevertex, baseinstance )
%
% C function:  void glDrawElementsInstancedBaseVertexBaseInstance(GLenum mode, GLsizei count, GLenum type, const void* indices, GLsizei instancecount, GLint basevertex, GLuint baseinstance)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glDrawElementsInstancedBaseVertexBaseInstance', mode, count, type, indices, instancecount, basevertex, baseinstance );

return
