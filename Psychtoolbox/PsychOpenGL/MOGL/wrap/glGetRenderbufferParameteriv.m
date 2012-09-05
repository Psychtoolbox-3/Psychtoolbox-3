function params = glGetRenderbufferParameteriv( target, pname )

% glGetRenderbufferParameteriv  Interface to OpenGL function glGetRenderbufferParameteriv
%
% usage:  params = glGetRenderbufferParameteriv( target, pname )
%
% C function:  void glGetRenderbufferParameteriv(GLenum target, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params=int32(repmat(intmax, [ 16 1 ]));
moglcore( 'glGetRenderbufferParameteriv', target, pname, params );
params = params(find(params~=intmax));

return
