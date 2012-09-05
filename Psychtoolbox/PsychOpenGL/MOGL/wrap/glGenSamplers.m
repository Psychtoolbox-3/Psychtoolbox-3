function samplers = glGenSamplers( count )

% glGenSamplers  Interface to OpenGL function glGenSamplers
%
% usage:  samplers = glGenSamplers( count )
%
% C function:  void glGenSamplers(GLsizei count, GLuint* samplers)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

samplers = uint32(zeros(1,count));

moglcore( 'glGenSamplers', count, samplers );

return
