function r = glIsProgramPipeline( pipeline )

% glIsProgramPipeline  Interface to OpenGL function glIsProgramPipeline
%
% usage:  r = glIsProgramPipeline( pipeline )
%
% C function:  GLboolean glIsProgramPipeline(GLuint pipeline)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

r = moglcore( 'glIsProgramPipeline', pipeline );

return
