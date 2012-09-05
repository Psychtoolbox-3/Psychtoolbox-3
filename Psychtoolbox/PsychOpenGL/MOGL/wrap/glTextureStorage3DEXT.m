function glTextureStorage3DEXT( texture, target, levels, internalformat, width, height, depth )

% glTextureStorage3DEXT  Interface to OpenGL function glTextureStorage3DEXT
%
% usage:  glTextureStorage3DEXT( texture, target, levels, internalformat, width, height, depth )
%
% C function:  void glTextureStorage3DEXT(GLuint texture, GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glTextureStorage3DEXT', texture, target, levels, internalformat, width, height, depth );

return
