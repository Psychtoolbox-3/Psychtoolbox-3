function glBlendEquationiARB( buf, mode )

% glBlendEquationiARB  Interface to OpenGL function glBlendEquationiARB
%
% usage:  glBlendEquationiARB( buf, mode )
%
% C function:  void glBlendEquationiARB(GLuint buf, GLenum mode)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glBlendEquationiARB', buf, mode );

return
