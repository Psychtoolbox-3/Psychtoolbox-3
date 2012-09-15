function [ length, size, type, name ] = glGetActiveAttrib( program, index, bufSize )

% glGetActiveAttrib  Interface to OpenGL function glGetActiveAttrib
%
% usage:  [ length, size, type, name ] = glGetActiveAttrib( program, index, bufSize )
%
% C function:  void glGetActiveAttrib(GLuint program, GLuint index, GLsizei bufSize, GLsizei* length, GLsizei* size, GLenum* type, GLchar* name)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

length = int32(0);
size = int32(0);
type = uint32(0);
name = uint8(zeros(1,bufSize));

moglcore( 'glGetActiveAttrib', program, index, bufSize, length, size, type, name );
name = char(name);

return
