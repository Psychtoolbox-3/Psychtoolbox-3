function [ length, label ] = glGetObjectPtrLabel( ptr, bufSize )

% glGetObjectPtrLabel  Interface to OpenGL function glGetObjectPtrLabel
%
% usage:  [ length, label ] = glGetObjectPtrLabel( ptr, bufSize )
%
% C function:  void glGetObjectPtrLabel(const void* ptr, GLsizei bufSize, GLsizei* length, GLchar* label)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

length = int32(0);
label = uint8(zeros(1, bufSize));

moglcore( 'glGetObjectPtrLabel', ptr, bufSize, length, label );
label = char(label);

return
