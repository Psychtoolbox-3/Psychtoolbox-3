function glDeleteRenderbuffers( n, renderbuffers )

% glDeleteRenderbuffers  Interface to OpenGL function glDeleteRenderbuffers
%
% usage:  glDeleteRenderbuffers( n, renderbuffers )
%
% C function:  void glDeleteRenderbuffers(GLsizei n, const GLuint* renderbuffers)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteRenderbuffers', n, uint32(renderbuffers) );

return
