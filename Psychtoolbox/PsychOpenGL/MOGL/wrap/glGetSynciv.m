function [ length, values ] = glGetSynciv( sync, pname, bufSize )

% glGetSynciv  Interface to OpenGL function glGetSynciv
%
% usage:  [ length, values ] = glGetSynciv( sync, pname, bufSize )
%
% C function:  void glGetSynciv(GLsync sync, GLenum pname, GLsizei bufSize, GLsizei* length, GLint* values)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

if ~strcmp(class(sync),'double'),
	error([ 'argument ''sync'' must be a pointer coded as type double ' ]);
end

length = int32(0);
values = int32(zeros(1, bufSize));

moglcore( 'glGetSynciv', sync, pname, bufSize, length, values );

return
