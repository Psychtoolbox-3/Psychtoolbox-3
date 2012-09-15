function glDrawArraysInstancedARB( mode, first, count, primcount )

% glDrawArraysInstancedARB  Interface to OpenGL function glDrawArraysInstancedARB
%
% usage:  glDrawArraysInstancedARB( mode, first, count, primcount )
%
% C function:  void glDrawArraysInstancedARB(GLenum mode, GLint first, GLsizei count, GLsizei primcount)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glDrawArraysInstancedARB', mode, first, count, primcount );

return
