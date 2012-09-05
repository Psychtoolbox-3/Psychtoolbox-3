function glObjectPtrLabel( ptr, length, label )

% glObjectPtrLabel  Interface to OpenGL function glObjectPtrLabel
%
% usage:  glObjectPtrLabel( ptr, length, label )
%
% C function:  void glObjectPtrLabel(const void* ptr, GLsizei length, const GLchar* label)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glObjectPtrLabel', ptr, length, uint8(label) );

return
