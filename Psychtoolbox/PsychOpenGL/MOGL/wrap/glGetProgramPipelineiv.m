function params = glGetProgramPipelineiv( pipeline, pname )

% glGetProgramPipelineiv  Interface to OpenGL function glGetProgramPipelineiv
%
% usage:  params = glGetProgramPipelineiv( pipeline, pname )
%
% C function:  void glGetProgramPipelineiv(GLuint pipeline, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params=int32(repmat(intmax, [ 16 1 ]));
moglcore( 'glGetProgramPipelineiv', pipeline, pname, params );
params = params(find(params~=intmax));

return
