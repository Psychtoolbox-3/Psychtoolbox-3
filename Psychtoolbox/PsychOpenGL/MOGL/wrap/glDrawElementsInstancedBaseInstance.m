function glDrawElementsInstancedBaseInstance( mode, count, type, indices, instancecount, baseinstance )

% glDrawElementsInstancedBaseInstance  Interface to OpenGL function glDrawElementsInstancedBaseInstance
%
% usage:  glDrawElementsInstancedBaseInstance( mode, count, type, indices, instancecount, baseinstance )
%
% C function:  void glDrawElementsInstancedBaseInstance(GLenum mode, GLsizei count, GLenum type, const void* indices, GLsizei instancecount, GLuint baseinstance)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glDrawElementsInstancedBaseInstance', mode, count, type, indices, instancecount, baseinstance );

return
