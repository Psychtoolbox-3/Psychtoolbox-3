function glBindProgramPipeline( pipeline )

% glBindProgramPipeline  Interface to OpenGL function glBindProgramPipeline
%
% usage:  glBindProgramPipeline( pipeline )
%
% C function:  void glBindProgramPipeline(GLuint pipeline)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glBindProgramPipeline', pipeline );

return
