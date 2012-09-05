function glDeleteFramebuffers( n, framebuffers )

% glDeleteFramebuffers  Interface to OpenGL function glDeleteFramebuffers
%
% usage:  glDeleteFramebuffers( n, framebuffers )
%
% C function:  void glDeleteFramebuffers(GLsizei n, const GLuint* framebuffers)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteFramebuffers', n, uint32(framebuffers) );

return
