function glTexStorage3DMultisample( target, samples, internalformat, width, height, depth, fixedsamplelocations )

% glTexStorage3DMultisample  Interface to OpenGL function glTexStorage3DMultisample
%
% usage:  glTexStorage3DMultisample( target, samples, internalformat, width, height, depth, fixedsamplelocations )
%
% C function:  void glTexStorage3DMultisample(GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glTexStorage3DMultisample', target, samples, internalformat, width, height, depth, fixedsamplelocations );

return
