function [ length, infoLog ] = glGetShaderInfoLog( shader, bufSize )

% glGetShaderInfoLog  Interface to OpenGL function glGetShaderInfoLog
%
% usage:  [ length, infoLog ] = glGetShaderInfoLog( shader, bufSize )
%
% C function:  void glGetShaderInfoLog(GLuint shader, GLsizei bufSize, GLsizei* length, GLchar* infoLog)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
infoLog = uint8(zeros(1,bufSize));

moglcore( 'glGetShaderInfoLog', shader, bufSize, length, infoLog );
infoLog = char(infoLog);

return
