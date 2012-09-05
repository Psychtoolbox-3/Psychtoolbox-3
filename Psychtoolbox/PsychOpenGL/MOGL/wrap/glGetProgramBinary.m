function [ length, binaryFormat, binary ] = glGetProgramBinary( program, bufSize )

% glGetProgramBinary  Interface to OpenGL function glGetProgramBinary
%
% usage:  [ length, binaryFormat, binary ] = glGetProgramBinary( program, bufSize )
%
% C function:  void glGetProgramBinary(GLuint program, GLsizei bufSize, GLsizei* length, GLenum* binaryFormat, GLvoid* binary)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
binaryFormat = uint32(0);
binary = uint8(zeros(1,bufSize));

moglcore( 'glGetProgramBinary', program, bufSize, length, binaryFormat, binary );

return
