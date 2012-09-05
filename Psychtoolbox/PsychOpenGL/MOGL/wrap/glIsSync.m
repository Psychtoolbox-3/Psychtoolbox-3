function r = glIsSync( sync )

% glIsSync  Interface to OpenGL function glIsSync
%
% usage:  r = glIsSync( sync )
%
% C function:  GLboolean glIsSync(GLsync sync)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

if ~strcmp(class(sync),'double'),
	error([ 'argument ''sync'' must be a pointer coded as type double ' ]);
end

r = moglcore( 'glIsSync', sync );

return
