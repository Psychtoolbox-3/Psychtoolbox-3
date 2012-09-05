function params = glGetSamplerParameterfv( sampler, pname )

% glGetSamplerParameterfv  Interface to OpenGL function glGetSamplerParameterfv
%
% usage:  params = glGetSamplerParameterfv( sampler, pname )
%
% C function:  void glGetSamplerParameterfv(GLuint sampler, GLenum pname, GLfloat* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = single(repmat(NaN,[ 16 1 ]));
moglcore( 'glGetSamplerParameterfv', sampler, pname, params );
params = double(params);
params = params(find(~isnan(params))); %#ok<FNDSB>

return
