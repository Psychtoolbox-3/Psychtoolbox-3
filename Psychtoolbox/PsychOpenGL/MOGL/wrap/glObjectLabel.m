function glObjectLabel( identifier, name, length, label )

% glObjectLabel  Interface to OpenGL function glObjectLabel
%
% usage:  glObjectLabel( identifier, name, length, label )
%
% C function:  void glObjectLabel(GLenum identifier, GLuint name, GLsizei length, const GLchar* label)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glObjectLabel', identifier, name, length, uint8(label) );

return
