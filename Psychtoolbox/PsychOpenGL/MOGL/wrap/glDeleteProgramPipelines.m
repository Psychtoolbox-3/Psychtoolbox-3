function glDeleteProgramPipelines( n, pipelines )

% glDeleteProgramPipelines  Interface to OpenGL function glDeleteProgramPipelines
%
% usage:  glDeleteProgramPipelines( n, pipelines )
%
% C function:  void glDeleteProgramPipelines(GLsizei n, const GLuint* pipelines)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteProgramPipelines', n, uint32(pipelines) );

return
