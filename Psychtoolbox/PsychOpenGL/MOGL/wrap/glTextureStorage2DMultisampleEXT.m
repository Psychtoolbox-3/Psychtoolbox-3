function glTextureStorage2DMultisampleEXT( texture, target, samples, internalformat, width, height, fixedsamplelocations )

% glTextureStorage2DMultisampleEXT  Interface to OpenGL function glTextureStorage2DMultisampleEXT
%
% usage:  glTextureStorage2DMultisampleEXT( texture, target, samples, internalformat, width, height, fixedsamplelocations )
%
% C function:  void glTextureStorage2DMultisampleEXT(GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glTextureStorage2DMultisampleEXT', texture, target, samples, internalformat, width, height, fixedsamplelocations );

return
