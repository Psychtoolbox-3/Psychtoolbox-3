function glMultiDrawElements( mode, count, type, indices, drawcount)

% glMultiDrawElements  Interface to OpenGL function glMultiDrawElements
%
% usage:  glMultiDrawElements( mode, count, type, indices, drawcount)
%
% Note: indices must be a cell array whose cells contain uint32() type
% vectors with the indices!
%
% C function:  void glMultiDrawElements(GLenum mode, const GLsizei* count, GLenum type, const GLvoid** indices, GLsizei drawcount)

% 1-Sep-2012 -- created (generated manually)

% ---protected---

if nargin~=5,
    error('invalid number of arguments');
end

% Ok, the moglcore() implementation doesn't work due to the double-indirect
% pointer, indices[] being an array of pointers to actual separate arrays
% of indices -- Difficult to handle / to pass from Matlab API to C API.
%
% Instead, we reimplement the behaviour of the function by iterative use of
% glDrawElements(), following the specification from the official
% spec of glMultiDrawElements(), which says exactly how to do
% this / defines the function in terms of iterating:
for i = 1:drawcount
    if count(i) > 0
        glDrawElements(mode, count(i), type, uint32(indices{i}));
    end
end

return
