function [ length, params ] = glGetProgramResourceiv( program, programInterface, index, propCount, props, bufSize )

% glGetProgramResourceiv  Interface to OpenGL function glGetProgramResourceiv
%
% usage:  [ length, params ] = glGetProgramResourceiv( program, programInterface, index, propCount, props, bufSize )
%
% C function:  void glGetProgramResourceiv(GLuint program, GLenum programInterface, GLuint index, GLsizei propCount, const GLenum* props, GLsizei bufSize, GLsizei* length, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=6,
    error('invalid number of arguments');
end

length = int32(0);
params = int32(zeros(1, bufSize);

moglcore( 'glGetProgramResourceiv', program, programInterface, index, propCount, uint32(props), bufSize, length, params );

return
