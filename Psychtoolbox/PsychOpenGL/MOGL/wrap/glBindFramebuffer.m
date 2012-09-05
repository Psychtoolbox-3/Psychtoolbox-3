function glBindFramebuffer( target, framebuffer )

% glBindFramebuffer  Interface to OpenGL function glBindFramebuffer
%
% usage:  glBindFramebuffer( target, framebuffer )
%
% C function:  void glBindFramebuffer(GLenum target, GLuint framebuffer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glBindFramebuffer', target, framebuffer );

return
