function params = glGetInternalformati64v( target, internalformat, pname, bufSize )

% glGetInternalformati64v  Interface to OpenGL function glGetInternalformati64v
%
% usage:  params = glGetInternalformati64v( target, internalformat, pname, bufSize )
%
% C function:  void glGetInternalformati64v(GLenum target, GLenum internalformat, GLenum pname, GLsizei bufSize, GLint64* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

params = int64(zeros(1,bufSize));

moglcore( 'glGetInternalformati64v', target, internalformat, pname, bufSize, params );

return
