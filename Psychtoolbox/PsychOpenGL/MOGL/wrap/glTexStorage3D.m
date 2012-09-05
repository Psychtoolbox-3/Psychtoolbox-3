function glTexStorage3D( target, levels, internalformat, width, height, depth )

% glTexStorage3D  Interface to OpenGL function glTexStorage3D
%
% usage:  glTexStorage3D( target, levels, internalformat, width, height, depth )
%
% C function:  void glTexStorage3D(GLenum target, GLsizei levels, GLenum internalformat, GLsizei width, GLsizei height, GLsizei depth)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=6,
    error('invalid number of arguments');
end

moglcore( 'glTexStorage3D', target, levels, internalformat, width, height, depth );

return
