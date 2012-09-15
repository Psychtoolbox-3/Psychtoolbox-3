function params = glGetSamplerParameterIuiv( sampler, pname )

% glGetSamplerParameterIuiv  Interface to OpenGL function glGetSamplerParameterIuiv
%
% usage:  params = glGetSamplerParameterIuiv( sampler, pname )
%
% C function:  void glGetSamplerParameterIuiv(GLuint sampler, GLenum pname, GLuint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params=uint32(repmat(intmax('uint32'), [ 16 1 ]));
moglcore( 'glGetSamplerParameterIuiv', sampler, pname, params );
params = params(find(params~=intmax('uint32')));
return
