function params = glGetNamedStringivARB( namelen, name, pname )

% glGetNamedStringivARB  Interface to OpenGL function glGetNamedStringivARB
%
% usage:  params = glGetNamedStringivARB( namelen, name, pname )
%
% C function:  void glGetNamedStringivARB(GLint namelen, const GLchar* name, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetNamedStringivARB', namelen, uint8(name), pname, params );

return
