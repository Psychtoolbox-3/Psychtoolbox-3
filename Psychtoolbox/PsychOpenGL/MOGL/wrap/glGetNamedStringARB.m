function [ stringlen, string ] = glGetNamedStringARB( namelen, name, bufSize )

% glGetNamedStringARB  Interface to OpenGL function glGetNamedStringARB
%
% usage:  [ stringlen, string ] = glGetNamedStringARB( namelen, name, bufSize )
%
% C function:  void glGetNamedStringARB(GLint namelen, const GLchar* name, GLsizei bufSize, GLint* stringlen, GLchar* string)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

stringlen = int32(0);
string = uint8(zeros(1,bufSize));

moglcore( 'glGetNamedStringARB', namelen, uint8(name), bufSize, stringlen, string );
string = char(string);

return
