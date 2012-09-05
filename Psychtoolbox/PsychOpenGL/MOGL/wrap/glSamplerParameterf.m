function glSamplerParameterf( sampler, pname, param )

% glSamplerParameterf  Interface to OpenGL function glSamplerParameterf
%
% usage:  glSamplerParameterf( sampler, pname, param )
%
% C function:  void glSamplerParameterf(GLuint sampler, GLenum pname, GLfloat param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glSamplerParameterf', sampler, pname, param );

return
