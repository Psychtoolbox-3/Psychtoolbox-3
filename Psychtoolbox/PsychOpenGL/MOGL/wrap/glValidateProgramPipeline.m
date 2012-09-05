function glValidateProgramPipeline( pipeline )

% glValidateProgramPipeline  Interface to OpenGL function glValidateProgramPipeline
%
% usage:  glValidateProgramPipeline( pipeline )
%
% C function:  void glValidateProgramPipeline(GLuint pipeline)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glValidateProgramPipeline', pipeline );

return
