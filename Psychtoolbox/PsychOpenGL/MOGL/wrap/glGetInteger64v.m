function params = glGetInteger64v( pname )

% glGetInteger64v  Interface to OpenGL function glGetInteger64v
%
% usage:  params = glGetInteger64v( pname )
%
% C function:  void glGetInteger64v(GLenum pname, GLint64* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=1,
    error('invalid number of arguments');
end

% Allocate enough bufferspace for the biggest returned vector:
params=int64(repmat(intmax('int64'), [ 16 1 ]));

moglcore( 'glGetInteger64v', pname, params );
params = params(find(params~=intmax('int64')));

return
