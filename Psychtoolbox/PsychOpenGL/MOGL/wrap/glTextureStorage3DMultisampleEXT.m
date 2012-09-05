function glTextureStorage3DMultisampleEXT( texture, target, samples, internalformat, width, height, depth, fixedsamplelocations )

% glTextureStorage3DMultisampleEXT  Interface to OpenGL function glTextureStorage3DMultisampleEXT
%
% usage:  glTextureStorage3DMultisampleEXT( texture, target, samples, internalformat, width, height, depth, fixedsamplelocations )
%
% C function:  void glTextureStorage3DMultisampleEXT(GLuint texture, GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=8,
    error('invalid number of arguments');
end

moglcore( 'glTextureStorage3DMultisampleEXT', texture, target, samples, internalformat, width, height, depth, fixedsamplelocations );

return
