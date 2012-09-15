function params = glGetQueryiv( target, pname )

% glGetQueryiv  Interface to OpenGL function glGetQueryiv
%
% usage:  params = glGetQueryiv( target, pname )
%
% C function:  void glGetQueryiv(GLenum target, GLenum pname, GLint* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetQueryiv', target, pname, params );

return
