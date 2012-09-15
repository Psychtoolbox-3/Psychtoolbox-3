function params = glGetNamedFramebufferParameterivEXT( framebuffer, pname )

% glGetNamedFramebufferParameterivEXT  Interface to OpenGL function glGetNamedFramebufferParameterivEXT
%
% usage:  params = glGetNamedFramebufferParameterivEXT( framebuffer, pname )
%
% C function:  void glGetNamedFramebufferParameterivEXT(GLuint framebuffer, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetNamedFramebufferParameterivEXT', framebuffer, pname, params );

return
