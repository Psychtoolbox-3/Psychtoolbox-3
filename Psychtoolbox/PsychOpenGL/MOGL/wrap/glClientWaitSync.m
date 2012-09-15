function r = glClientWaitSync( sync, flags, timeout )

% glClientWaitSync  Interface to OpenGL function glClientWaitSync
%
% usage:  r = glClientWaitSync( sync, flags, timeout )
%
% C function:  GLenum glClientWaitSync(GLsync sync, GLbitfield flags, GLuint64 timeout)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

if ~strcmp(class(sync),'double'),
	error([ 'argument ''sync'' must be a pointer coded as type double ' ]);
end

r = moglcore( 'glClientWaitSync', sync, flags, uint64(timeout) );

return
