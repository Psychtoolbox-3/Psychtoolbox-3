function params = glGetUniformSubroutineuiv( shadertype, location )

% glGetUniformSubroutineuiv  Interface to OpenGL function glGetUniformSubroutineuiv
%
% usage:  params = glGetUniformSubroutineuiv( shadertype, location )
%
% C function:  void glGetUniformSubroutineuiv(GLenum shadertype, GLint location, GLuint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

% Allocate enough bufferspace for the biggest returned vector:
params=uint32(repmat(intmax('uint32'), [ 16 1 ]));

moglcore( 'glGetUniformSubroutineuiv', shadertype, location, params );
params = params(find(params~=intmax('uint32')));

return
