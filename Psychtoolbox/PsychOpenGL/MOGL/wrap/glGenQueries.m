function ids = glGenQueries( n )

% glGenQueries  Interface to OpenGL function glGenQueries
%
% usage:  ids = glGenQueries( n )
%
% C function:  void glGenQueries(GLsizei n, GLuint* ids)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

ids = uint32(zeros(1,n));

moglcore( 'glGenQueries', n, ids );

return
