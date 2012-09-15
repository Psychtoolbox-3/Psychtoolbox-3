function glDrawRangeElementsBaseVertex( mode, start, end, count, type, indices, basevertex )

% glDrawRangeElementsBaseVertex  Interface to OpenGL function glDrawRangeElementsBaseVertex
%
% usage:  glDrawRangeElementsBaseVertex( mode, start, end, count, type, indices, basevertex )
%
% C function:  void glDrawRangeElementsBaseVertex(GLenum mode, GLuint start, GLuint end, GLsizei count, GLenum type, const GLvoid* indices, GLint basevertex)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glDrawRangeElementsBaseVertex', mode, start, end, count, type, indices, basevertex );

return
