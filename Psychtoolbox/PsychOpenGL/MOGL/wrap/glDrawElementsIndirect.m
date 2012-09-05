function glDrawElementsIndirect( mode, type, indirect )

% glDrawElementsIndirect  Interface to OpenGL function glDrawElementsIndirect
%
% usage:  glDrawElementsIndirect( mode, type, indirect )
%
% C function:  void glDrawElementsIndirect(GLenum mode, GLenum type, const GLvoid* indirect)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDrawElementsIndirect', mode, type, indirect );

return
