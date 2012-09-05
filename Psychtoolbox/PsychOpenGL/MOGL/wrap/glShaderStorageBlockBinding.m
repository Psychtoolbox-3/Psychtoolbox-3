function glShaderStorageBlockBinding( program, storageBlockIndex, storageBlockBinding )

% glShaderStorageBlockBinding  Interface to OpenGL function glShaderStorageBlockBinding
%
% usage:  glShaderStorageBlockBinding( program, storageBlockIndex, storageBlockBinding )
%
% C function:  void glShaderStorageBlockBinding(GLuint program, GLuint storageBlockIndex, GLuint storageBlockBinding)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glShaderStorageBlockBinding', program, storageBlockIndex, storageBlockBinding );

return
