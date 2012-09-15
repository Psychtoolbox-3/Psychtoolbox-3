function glScissorArrayv( first, count, v )

% glScissorArrayv  Interface to OpenGL function glScissorArrayv
%
% usage:  glScissorArrayv( first, count, v )
%
% C function:  void glScissorArrayv(GLuint first, GLsizei count, const GLint* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glScissorArrayv', first, count, int32(v) );

return
