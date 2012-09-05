function renderbuffers = glGenRenderbuffers( n )

% glGenRenderbuffers  Interface to OpenGL function glGenRenderbuffers
%
% usage:  renderbuffers = glGenRenderbuffers( n )
%
% C function:  void glGenRenderbuffers(GLsizei n, GLuint* renderbuffers)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

renderbuffers = uint32(zeros(1,n));

moglcore( 'glGenRenderbuffers', n, renderbuffers );

return
