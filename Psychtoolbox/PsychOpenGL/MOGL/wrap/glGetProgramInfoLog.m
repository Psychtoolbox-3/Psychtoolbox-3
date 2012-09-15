function [ length, infoLog ] = glGetProgramInfoLog( program, bufSize )

% glGetProgramInfoLog  Interface to OpenGL function glGetProgramInfoLog
%
% usage:  [ length, infoLog ] = glGetProgramInfoLog( program, bufSize )
%
% C function:  void glGetProgramInfoLog(GLuint program, GLsizei bufSize, GLsizei* length, GLchar* infoLog)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
infoLog = uint8(zeros(1,bufSize));

moglcore( 'glGetProgramInfoLog', program, bufSize, length, infoLog );
infoLog = char(infoLog);

return
