function data = glGetDoublei_v( target, index )

% glGetDoublei_v  Interface to OpenGL function glGetDoublei_v
%
% usage:  data = glGetDoublei_v( target, index )
%
% C function:  void glGetDoublei_v(GLenum target, GLuint index, GLdouble* data)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

data = double(0);

moglcore( 'glGetDoublei_v', target, index, data );

return
