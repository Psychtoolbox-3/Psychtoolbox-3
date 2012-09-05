function glFramebufferTexture3D( target, attachment, textarget, texture, level, zoffset )

% glFramebufferTexture3D  Interface to OpenGL function glFramebufferTexture3D
%
% usage:  glFramebufferTexture3D( target, attachment, textarget, texture, level, zoffset )
%
% C function:  void glFramebufferTexture3D(GLenum target, GLenum attachment, GLenum textarget, GLuint texture, GLint level, GLint zoffset)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferTexture3D', target, attachment, textarget, texture, level, zoffset );

return
