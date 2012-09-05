function glMultiDrawArraysIndirect( mode, indirect, drawcount, stride )

% glMultiDrawArraysIndirect  Interface to OpenGL function glMultiDrawArraysIndirect
%
% usage:  glMultiDrawArraysIndirect( mode, indirect, drawcount, stride )
%
% C function:  void glMultiDrawArraysIndirect(GLenum mode, const void* indirect, GLsizei drawcount, GLsizei stride)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glMultiDrawArraysIndirect', mode, indirect, drawcount, stride );

return
