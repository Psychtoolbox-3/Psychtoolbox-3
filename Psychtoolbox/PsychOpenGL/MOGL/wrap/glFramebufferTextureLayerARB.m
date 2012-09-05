function glFramebufferTextureLayerARB( target, attachment, texture, level, layer )

% glFramebufferTextureLayerARB  Interface to OpenGL function glFramebufferTextureLayerARB
%
% usage:  glFramebufferTextureLayerARB( target, attachment, texture, level, layer )
%
% C function:  void glFramebufferTextureLayerARB(GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTextureLayerARB', target, attachment, texture, level, layer );

return
