function glEndQueryIndexed( target, index )

% glEndQueryIndexed  Interface to OpenGL function glEndQueryIndexed
%
% usage:  glEndQueryIndexed( target, index )
%
% C function:  void glEndQueryIndexed(GLenum target, GLuint index)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glEndQueryIndexed', target, index );

return
