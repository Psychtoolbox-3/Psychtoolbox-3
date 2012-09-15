function glTexImage3DMultisample( target, samples, internalformat, width, height, depth, fixedsamplelocations )

% glTexImage3DMultisample  Interface to OpenGL function glTexImage3DMultisample
%
% usage:  glTexImage3DMultisample( target, samples, internalformat, width, height, depth, fixedsamplelocations )
%
% C function:  void glTexImage3DMultisample(GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glTexImage3DMultisample', target, samples, internalformat, width, height, depth, fixedsamplelocations );

return
