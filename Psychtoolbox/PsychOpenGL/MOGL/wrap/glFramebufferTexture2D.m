function glFramebufferTexture2D( target, attachment, textarget, texture, level )

% glFramebufferTexture2D  Interface to OpenGL function glFramebufferTexture2D
%
% usage:  glFramebufferTexture2D( target, attachment, textarget, texture, level )
%
% C function:  void glFramebufferTexture2D(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTexture2D', target, attachment, textarget, texture, level );

return
