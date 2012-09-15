function glTexStorage2DMultisample( target, samples, internalformat, width, height, fixedsamplelocations )

% glTexStorage2DMultisample  Interface to OpenGL function glTexStorage2DMultisample
%
% usage:  glTexStorage2DMultisample( target, samples, internalformat, width, height, fixedsamplelocations )
%
% C function:  void glTexStorage2DMultisample(GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glTexStorage2DMultisample', target, samples, internalformat, width, height, fixedsamplelocations );

return
