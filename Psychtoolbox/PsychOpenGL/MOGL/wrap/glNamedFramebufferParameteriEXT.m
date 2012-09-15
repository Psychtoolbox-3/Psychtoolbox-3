function glNamedFramebufferParameteriEXT( framebuffer, pname, param )

% glNamedFramebufferParameteriEXT  Interface to OpenGL function glNamedFramebufferParameteriEXT
%
% usage:  glNamedFramebufferParameteriEXT( framebuffer, pname, param )
%
% C function:  void glNamedFramebufferParameteriEXT(GLuint framebuffer, GLenum pname, GLint param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glNamedFramebufferParameteriEXT', framebuffer, pname, param );

return
