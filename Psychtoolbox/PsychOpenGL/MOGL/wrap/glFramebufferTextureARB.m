function glFramebufferTextureARB( target, attachment, texture, level )

% glFramebufferTextureARB  Interface to OpenGL function glFramebufferTextureARB
%
% usage:  glFramebufferTextureARB( target, attachment, texture, level )
%
% C function:  void glFramebufferTextureARB(GLenum target, GLenum attachment, GLuint texture, GLint level)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTextureARB', target, attachment, texture, level );

return
