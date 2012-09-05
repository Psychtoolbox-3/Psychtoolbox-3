function glGenerateMipmap( target )

% glGenerateMipmap  Interface to OpenGL function glGenerateMipmap
%
% usage:  glGenerateMipmap( target )
%
% C function:  void glGenerateMipmap(GLenum target)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glGenerateMipmap', target );

return
