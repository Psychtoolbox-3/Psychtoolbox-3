function ids = glGenTransformFeedbacks( n )

% glGenTransformFeedbacks  Interface to OpenGL function glGenTransformFeedbacks
%
% usage:  ids = glGenTransformFeedbacks( n )
%
% C function:  void glGenTransformFeedbacks(GLsizei n, GLuint* ids)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

ids = uint32(zeros(1,n));

moglcore( 'glGenTransformFeedbacks', n, ids );

return
