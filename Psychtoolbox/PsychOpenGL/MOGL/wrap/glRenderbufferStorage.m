function glRenderbufferStorage( target, internalformat, width, height )

% glRenderbufferStorage  Interface to OpenGL function glRenderbufferStorage
%
% usage:  glRenderbufferStorage( target, internalformat, width, height )
%
% C function:  void glRenderbufferStorage(GLenum target, GLenum internalformat, GLsizei width, GLsizei height)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glRenderbufferStorage', target, internalformat, width, height );

return
