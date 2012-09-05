function glShaderBinary( count, shaders, binaryformat, binary, length )

% glShaderBinary  Interface to OpenGL function glShaderBinary
%
% usage:  glShaderBinary( count, shaders, binaryformat, binary, length )
%
% C function:  void glShaderBinary(GLsizei count, const GLuint* shaders, GLenum binaryformat, const GLvoid* binary, GLsizei length)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glShaderBinary', count, uint32(shaders), binaryformat, binary, length );

return
