function glPatchParameteri( pname, value )

% glPatchParameteri  Interface to OpenGL function glPatchParameteri
%
% usage:  glPatchParameteri( pname, value )
%
% C function:  void glPatchParameteri(GLenum pname, GLint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glPatchParameteri', pname, value );

return
