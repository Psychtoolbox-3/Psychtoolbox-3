function glMultiDrawElementsBaseVertex( mode, count, type, indices, drawcount, basevertex )

% glMultiDrawElementsBaseVertex  Interface to OpenGL function glMultiDrawElementsBaseVertex
%
% usage:  glMultiDrawElementsBaseVertex( mode, count, type, indices, drawcount, basevertex )
%
% Note: indices must be a cell array whose cells contain uint32() type
% vectors with the indices!
%
% C function:  void glMultiDrawElementsBaseVertex(GLenum mode, const GLsizei* count, GLenum type, const GLvoid** indices, GLsizei drawcount, const GLint* basevertex)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---protected---

if nargin~=6,
    error('invalid number of arguments');
end

% Ok, the moglcore() implementation doesn't work due to the double-indirect
% pointer, indices[] being an array of pointers to actual separate arrays
% of indices -- Difficult to handle / to pass from Matlab API to C API.
%
% Therefore this call is disabled:
% moglcore( 'glMultiDrawElementsBaseVertex', mode, int32(count), type, const, drawcount, int32(basevertex) );

% Instead, we reimplement the behaviour of the function by iterative use of
% glDrawElementsBaseVertex(), following the specification from the official
% spec of glMultiDrawElementsBaseVertex(), which says exactly how to do
% this / defines the function in terms of iterating:
for i = 1:drawcount
    if count(i) > 0
        glDrawElementsBaseVertex(mode, count(i), type, uint32(indices{i}), basevertex(i));
    end
end

return
