function glDrawElementsBaseVertex( mode, count, type, indices, basevertex )

% glDrawElementsBaseVertex  Interface to OpenGL function glDrawElementsBaseVertex
%
% usage:  glDrawElementsBaseVertex( mode, count, type, indices, basevertex )
%
% C function:  void glDrawElementsBaseVertex(GLenum mode, GLsizei count, GLenum type, const GLvoid* indices, GLint basevertex)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glDrawElementsBaseVertex', mode, count, type, indices, basevertex );

return
