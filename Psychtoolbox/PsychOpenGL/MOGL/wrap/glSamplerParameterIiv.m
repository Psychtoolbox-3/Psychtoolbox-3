function glSamplerParameterIiv( sampler, pname, param )

% glSamplerParameterIiv  Interface to OpenGL function glSamplerParameterIiv
%
% usage:  glSamplerParameterIiv( sampler, pname, param )
%
% C function:  void glSamplerParameterIiv(GLuint sampler, GLenum pname, const GLint* param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glSamplerParameterIiv', sampler, pname, int32(param) );

return
