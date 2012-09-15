function glSamplerParameteriv( sampler, pname, param )

% glSamplerParameteriv  Interface to OpenGL function glSamplerParameteriv
%
% usage:  glSamplerParameteriv( sampler, pname, param )
%
% C function:  void glSamplerParameteriv(GLuint sampler, GLenum pname, const GLint* param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glSamplerParameteriv', sampler, pname, int32(param) );

return
