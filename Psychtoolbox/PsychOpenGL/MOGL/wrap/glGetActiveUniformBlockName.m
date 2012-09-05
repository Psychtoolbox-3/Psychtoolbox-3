function [ length, uniformBlockName ] = glGetActiveUniformBlockName( program, uniformBlockIndex, bufSize )

% glGetActiveUniformBlockName  Interface to OpenGL function glGetActiveUniformBlockName
%
% usage:  [ length, uniformBlockName ] = glGetActiveUniformBlockName( program, uniformBlockIndex, bufSize )
%
% C function:  void glGetActiveUniformBlockName(GLuint program, GLuint uniformBlockIndex, GLsizei bufSize, GLsizei* length, GLchar* uniformBlockName)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

length = int32(0);
uniformBlockName = uint8(zeros(1,bufSize));

moglcore( 'glGetActiveUniformBlockName', program, uniformBlockIndex, bufSize, length, uniformBlockName );
uniformBlockName = char(uniformBlockName);

return
