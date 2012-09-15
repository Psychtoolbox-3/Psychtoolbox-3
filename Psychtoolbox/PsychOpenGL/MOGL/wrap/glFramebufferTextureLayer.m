function glFramebufferTextureLayer( target, attachment, texture, level, layer )

% glFramebufferTextureLayer  Interface to OpenGL function glFramebufferTextureLayer
%
% usage:  glFramebufferTextureLayer( target, attachment, texture, level, layer )
%
% C function:  void glFramebufferTextureLayer(GLenum target, GLenum attachment, GLuint texture, GLint level, GLint layer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTextureLayer', target, attachment, texture, level, layer );

return
