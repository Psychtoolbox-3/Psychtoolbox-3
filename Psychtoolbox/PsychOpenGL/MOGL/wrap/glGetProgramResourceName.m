function [ length, name ] = glGetProgramResourceName( program, programInterface, index, bufSize )

% glGetProgramResourceName  Interface to OpenGL function glGetProgramResourceName
%
% usage:  [ length, name ] = glGetProgramResourceName( program, programInterface, index, bufSize )
%
% C function:  void glGetProgramResourceName(GLuint program, GLenum programInterface, GLuint index, GLsizei bufSize, GLsizei* length, GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

length = int32(0);
name = uint8(zeros(1,bufSize));

moglcore( 'glGetProgramResourceName', program, programInterface, index, bufSize, length, name );

return
