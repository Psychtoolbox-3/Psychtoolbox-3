function glFramebufferParameteri( target, pname, param )

% glFramebufferParameteri  Interface to OpenGL function glFramebufferParameteri
%
% usage:  glFramebufferParameteri( target, pname, param )
%
% C function:  void glFramebufferParameteri(GLenum target, GLenum pname, GLint param)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glFramebufferParameteri', target, pname, param );

return
