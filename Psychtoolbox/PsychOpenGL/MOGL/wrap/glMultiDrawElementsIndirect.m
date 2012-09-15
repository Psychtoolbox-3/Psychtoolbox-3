function glMultiDrawElementsIndirect( mode, type, indirect, drawcount, stride )

% glMultiDrawElementsIndirect  Interface to OpenGL function glMultiDrawElementsIndirect
%
% usage:  glMultiDrawElementsIndirect( mode, type, indirect, drawcount, stride )
%
% C function:  void glMultiDrawElementsIndirect(GLenum mode, GLenum type, const void* indirect, GLsizei drawcount, GLsizei stride)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glMultiDrawElementsIndirect', mode, type, indirect, drawcount, stride );

return
