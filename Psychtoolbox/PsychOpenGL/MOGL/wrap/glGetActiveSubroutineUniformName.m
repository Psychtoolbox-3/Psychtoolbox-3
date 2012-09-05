function [ length, name ] = glGetActiveSubroutineUniformName( program, shadertype, index, bufsize )

% glGetActiveSubroutineUniformName  Interface to OpenGL function glGetActiveSubroutineUniformName
%
% usage:  [ length, name ] = glGetActiveSubroutineUniformName( program, shadertype, index, bufsize )
%
% C function:  void glGetActiveSubroutineUniformName(GLuint program, GLenum shadertype, GLuint index, GLsizei bufsize, GLsizei* length, GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

length = int32(0);
name = uint8(zeros(1,bufsize));

moglcore( 'glGetActiveSubroutineUniformName', program, shadertype, index, bufsize, length, name );
name = char(name);

return
