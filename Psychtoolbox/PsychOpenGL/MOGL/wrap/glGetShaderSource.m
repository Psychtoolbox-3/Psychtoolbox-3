function [ length, source ] = glGetShaderSource( shader, bufSize )

% glGetShaderSource  Interface to OpenGL function glGetShaderSource
%
% usage:  [ length, source ] = glGetShaderSource( shader, bufSize )
%
% C function:  void glGetShaderSource(GLuint shader, GLsizei bufSize, GLsizei* length, GLchar* source)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
source = uint8(zeros(1,bufSize));

moglcore( 'glGetShaderSource', shader, bufSize, length, source );
source = char(source);

return
