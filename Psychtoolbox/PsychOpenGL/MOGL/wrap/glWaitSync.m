function glWaitSync( sync, flags, timeout )

% glWaitSync  Interface to OpenGL function glWaitSync
%
% usage:  glWaitSync( sync, flags, timeout )
%
% C function:  void glWaitSync(GLsync sync, GLbitfield flags, GLuint64 timeout)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

if ~strcmp(class(sync),'double'),
	error([ 'argument ''sync'' must be a pointer coded as type double ' ]);
end

moglcore( 'glWaitSync', sync, flags, uint64(timeout) );

return
