function glTexStorage1D( target, levels, internalformat, width )

% glTexStorage1D  Interface to OpenGL function glTexStorage1D
%
% usage:  glTexStorage1D( target, levels, internalformat, width )
%
% C function:  void glTexStorage1D(GLenum target, GLsizei levels, GLenum internalformat, GLsizei width)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glTexStorage1D', target, levels, internalformat, width );

return
