function glDrawElementsInstancedBaseVertex( mode, count, type, indices, instancecount, basevertex )

% glDrawElementsInstancedBaseVertex  Interface to OpenGL function glDrawElementsInstancedBaseVertex
%
% usage:  glDrawElementsInstancedBaseVertex( mode, count, type, indices, instancecount, basevertex )
%
% C function:  void glDrawElementsInstancedBaseVertex(GLenum mode, GLsizei count, GLenum type, const GLvoid* indices, GLsizei instancecount, GLint basevertex)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glDrawElementsInstancedBaseVertex', mode, count, type, indices, instancecount, basevertex );

return
