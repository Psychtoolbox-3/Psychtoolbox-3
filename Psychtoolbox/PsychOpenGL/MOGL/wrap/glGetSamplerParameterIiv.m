function params = glGetSamplerParameterIiv( sampler, pname )

% glGetSamplerParameterIiv  Interface to OpenGL function glGetSamplerParameterIiv
%
% usage:  params = glGetSamplerParameterIiv( sampler, pname )
%
% C function:  void glGetSamplerParameterIiv(GLuint sampler, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params=int32(repmat(intmax('int32'), [ 16 1 ]));
moglcore( 'glGetSamplerParameterIiv', sampler, pname, params );
params = params(find(params~=intmax('int32')));

return
