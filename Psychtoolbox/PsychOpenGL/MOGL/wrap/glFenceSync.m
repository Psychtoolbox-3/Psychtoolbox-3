function r = glFenceSync( condition, flags )

% glFenceSync  Interface to OpenGL function glFenceSync
%
% usage:  r = glFenceSync( condition, flags )
%
% C function:  GLsync glFenceSync(GLenum condition, GLbitfield flags)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

r = moglcore( 'glFenceSync', condition, flags );

return
