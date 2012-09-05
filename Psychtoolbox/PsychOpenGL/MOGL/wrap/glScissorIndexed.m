function glScissorIndexed( index, left, bottom, width, height )

% glScissorIndexed  Interface to OpenGL function glScissorIndexed
%
% usage:  glScissorIndexed( index, left, bottom, width, height )
%
% C function:  void glScissorIndexed(GLuint index, GLint left, GLint bottom, GLsizei width, GLsizei height)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glScissorIndexed', index, left, bottom, width, height );

return
