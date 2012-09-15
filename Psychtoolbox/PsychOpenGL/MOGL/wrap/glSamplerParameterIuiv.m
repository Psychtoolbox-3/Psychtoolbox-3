function glSamplerParameterIuiv( sampler, pname, param )

% glSamplerParameterIuiv  Interface to OpenGL function glSamplerParameterIuiv
%
% usage:  glSamplerParameterIuiv( sampler, pname, param )
%
% C function:  void glSamplerParameterIuiv(GLuint sampler, GLenum pname, const GLuint* param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glSamplerParameterIuiv', sampler, pname, uint32(param) );

return
