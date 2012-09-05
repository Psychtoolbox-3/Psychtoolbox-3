function glMemoryBarrier( barriers )

% glMemoryBarrier  Interface to OpenGL function glMemoryBarrier
%
% usage:  glMemoryBarrier( barriers )
%
% C function:  void glMemoryBarrier(GLbitfield barriers)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glMemoryBarrier', barriers );

return
