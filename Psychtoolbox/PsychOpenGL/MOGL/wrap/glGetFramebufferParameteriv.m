function params = glGetFramebufferParameteriv( target, pname )

% glGetFramebufferParameteriv  Interface to OpenGL function glGetFramebufferParameteriv
%
% usage:  params = glGetFramebufferParameteriv( target, pname )
%
% C function:  void glGetFramebufferParameteriv(GLenum target, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetFramebufferParameteriv', target, pname, params );

return
