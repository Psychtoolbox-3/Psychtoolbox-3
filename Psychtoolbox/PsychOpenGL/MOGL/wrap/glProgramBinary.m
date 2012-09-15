function glProgramBinary( program, binaryFormat, binary, length )

% glProgramBinary  Interface to OpenGL function glProgramBinary
%
% usage:  glProgramBinary( program, binaryFormat, binary, length )
%
% C function:  void glProgramBinary(GLuint program, GLenum binaryFormat, const GLvoid* binary, GLsizei length)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glProgramBinary', program, binaryFormat, binary, length );

return
