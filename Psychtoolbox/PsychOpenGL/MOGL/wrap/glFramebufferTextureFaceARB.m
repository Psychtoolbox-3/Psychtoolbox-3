function glFramebufferTextureFaceARB( target, attachment, texture, level, face )

% glFramebufferTextureFaceARB  Interface to OpenGL function glFramebufferTextureFaceARB
%
% usage:  glFramebufferTextureFaceARB( target, attachment, texture, level, face )
%
% C function:  void glFramebufferTextureFaceARB(GLenum target, GLenum attachment, GLuint texture, GLint level, GLenum face)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTextureFaceARB', target, attachment, texture, level, face );

return
