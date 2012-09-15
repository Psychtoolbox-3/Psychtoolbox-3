function params = glGetSamplerParameteriv( sampler, pname )

% glGetSamplerParameteriv  Interface to OpenGL function glGetSamplerParameteriv
%
% usage:  params = glGetSamplerParameteriv( sampler, pname )
%
% C function:  void glGetSamplerParameteriv(GLuint sampler, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params=int32(repmat(intmax, [ 16 1 ]));
moglcore( 'glGetSamplerParameteriv', sampler, pname, params );
params = params(find(params~=intmax));

return
