function glClearDepthf( d )

% glClearDepthf  Interface to OpenGL function glClearDepthf
%
% usage:  glClearDepthf( d )
%
% C function:  void glClearDepthf(GLfloat d)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glClearDepthf', d );

return
