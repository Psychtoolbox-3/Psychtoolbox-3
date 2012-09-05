function glTexImage2DMultisample( target, samples, internalformat, width, height, fixedsamplelocations )

% glTexImage2DMultisample  Interface to OpenGL function glTexImage2DMultisample
%
% usage:  glTexImage2DMultisample( target, samples, internalformat, width, height, fixedsamplelocations )
%
% C function:  void glTexImage2DMultisample(GLenum target, GLsizei samples, GLint internalformat, GLsizei width, GLsizei height, GLboolean fixedsamplelocations)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glTexImage2DMultisample', target, samples, internalformat, width, height, fixedsamplelocations );

return
