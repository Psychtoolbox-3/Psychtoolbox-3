function [ length, uniformName ] = glGetActiveUniformName( program, uniformIndex, bufSize )

% glGetActiveUniformName  Interface to OpenGL function glGetActiveUniformName
%
% usage:  [ length, uniformName ] = glGetActiveUniformName( program, uniformIndex, bufSize )
%
% C function:  void glGetActiveUniformName(GLuint program, GLuint uniformIndex, GLsizei bufSize, GLsizei* length, GLchar* uniformName)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

length = int32(0);
uniformName = uint8(zeros(1,bufSize));

moglcore( 'glGetActiveUniformName', program, uniformIndex, bufSize, length, uniformName );
uniformName = char(uniformName);

return
