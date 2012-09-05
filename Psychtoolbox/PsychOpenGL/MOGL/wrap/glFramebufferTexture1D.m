function glFramebufferTexture1D( target, attachment, textarget, texture, level )

% glFramebufferTexture1D  Interface to OpenGL function glFramebufferTexture1D
%
% usage:  glFramebufferTexture1D( target, attachment, textarget, texture, level )
%
% C function:  void glFramebufferTexture1D(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTexture1D', target, attachment, textarget, texture, level );

return
