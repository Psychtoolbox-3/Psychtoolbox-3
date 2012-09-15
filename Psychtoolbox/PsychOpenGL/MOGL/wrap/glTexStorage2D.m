function glTexStorage2D( target, levels, internalformat, width, height )

% glTexStorage2D  Interface to OpenGL function glTexStorage2D
%
% usage:  glTexStorage2D( target, levels, internalformat, width, height )
%
% C function:  void glTexStorage2D(GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glTexStorage2D', target, levels, internalformat, width, height );

return
