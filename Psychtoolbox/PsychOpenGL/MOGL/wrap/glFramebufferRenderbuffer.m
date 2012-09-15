function glFramebufferRenderbuffer( target, attachment, renderbuffertarget, renderbuffer )

% glFramebufferRenderbuffer  Interface to OpenGL function glFramebufferRenderbuffer
%
% usage:  glFramebufferRenderbuffer( target, attachment, renderbuffertarget, renderbuffer )
%
% C function:  void glFramebufferRenderbuffer(GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferRenderbuffer', target, attachment, renderbuffertarget, renderbuffer );

return
