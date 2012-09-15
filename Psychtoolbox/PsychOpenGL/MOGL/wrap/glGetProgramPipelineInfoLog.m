function [ length, infoLog ] = glGetProgramPipelineInfoLog( pipeline, bufSize )

% glGetProgramPipelineInfoLog  Interface to OpenGL function glGetProgramPipelineInfoLog
%
% usage:  [ length, infoLog ] = glGetProgramPipelineInfoLog( pipeline, bufSize )
%
% C function:  void glGetProgramPipelineInfoLog(GLuint pipeline, GLsizei bufSize, GLsizei* length, GLchar* infoLog)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
infoLog = uint8(zeros(1,bufSize));

moglcore( 'glGetProgramPipelineInfoLog', pipeline, bufSize, length, infoLog );
infoLog = char(infoLog);

return
