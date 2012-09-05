function params = glGetInternalformativ( target, internalformat, pname, bufSize )

% glGetInternalformativ  Interface to OpenGL function glGetInternalformativ
%
% usage:  params = glGetInternalformativ( target, internalformat, pname, bufSize )
%
% C function:  void glGetInternalformativ(GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

params = int32(zeros(1,bufSize));

moglcore( 'glGetInternalformativ', target, internalformat, pname, bufSize, params );

return
