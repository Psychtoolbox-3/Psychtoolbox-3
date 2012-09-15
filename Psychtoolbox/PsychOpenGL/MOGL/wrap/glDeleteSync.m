function glDeleteSync( sync )

% glDeleteSync  Interface to OpenGL function glDeleteSync
%
% usage:  glDeleteSync( sync )
%
% C function:  void glDeleteSync(GLsync sync)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

if ~strcmp(class(sync),'double'),
	error([ 'argument ''sync'' must be a pointer coded as type double ' ]);
end

moglcore( 'glDeleteSync', sync );

return
