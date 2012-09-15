function glTexBufferARB( target, internalformat, buffer )

% glTexBufferARB  Interface to OpenGL function glTexBufferARB
%
% usage:  glTexBufferARB( target, internalformat, buffer )
%
% C function:  void glTexBufferARB(GLenum target, GLenum internalformat, GLuint buffer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glTexBufferARB', target, internalformat, buffer );

return
