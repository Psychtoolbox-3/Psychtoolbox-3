function pipelines = glGenProgramPipelines( n )

% glGenProgramPipelines  Interface to OpenGL function glGenProgramPipelines
%
% usage:  pipelines = glGenProgramPipelines( n )
%
% C function:  void glGenProgramPipelines(GLsizei n, GLuint* pipelines)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

pipelines = uint32(zeros(1,n));

moglcore( 'glGenProgramPipelines', n, pipelines );

return
