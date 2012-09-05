function glRenderbufferStorageMultisample( target, samples, internalformat, width, height )

% glRenderbufferStorageMultisample  Interface to OpenGL function glRenderbufferStorageMultisample
%
% usage:  glRenderbufferStorageMultisample( target, samples, internalformat, width, height )
%
% C function:  void glRenderbufferStorageMultisample(GLenum target, GLsizei samples, GLenum internalformat, GLsizei width, GLsizei height)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glRenderbufferStorageMultisample', target, samples, internalformat, width, height );

return
