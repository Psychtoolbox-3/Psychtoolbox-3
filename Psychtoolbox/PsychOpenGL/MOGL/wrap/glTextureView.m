function glTextureView( texture, target, origtexture, internalformat, minlevel, numlevels, minlayer, numlayers )

% glTextureView  Interface to OpenGL function glTextureView
%
% usage:  glTextureView( texture, target, origtexture, internalformat, minlevel, numlevels, minlayer, numlayers )
%
% C function:  void glTextureView(GLuint texture, GLenum target, GLuint origtexture, GLenum internalformat, GLuint minlevel, GLuint numlevels, GLuint minlayer, GLuint numlayers)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=8,
    error('invalid number of arguments');
end

moglcore( 'glTextureView', texture, target, origtexture, internalformat, minlevel, numlevels, minlayer, numlayers );

return
