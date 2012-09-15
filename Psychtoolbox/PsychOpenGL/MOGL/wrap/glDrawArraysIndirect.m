function glDrawArraysIndirect( mode, indirect )

% glDrawArraysIndirect  Interface to OpenGL function glDrawArraysIndirect
%
% usage:  glDrawArraysIndirect( mode, indirect )
%
% C function:  void glDrawArraysIndirect(GLenum mode, const GLvoid* indirect)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDrawArraysIndirect', mode, indirect );

return
