function uniformIndices = glGetUniformIndices( program, uniformCount, const )

% glGetUniformIndices  Interface to OpenGL function glGetUniformIndices
%
% usage:  uniformIndices = glGetUniformIndices( program, uniformCount, const )
%
% C function:  void glGetUniformIndices(GLuint program, GLsizei uniformCount, const GLchar* const, GLuint* uniformIndices)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---

if nargin~=3,
    error('invalid number of arguments');
end

uniformIndices = uint32(0);

moglcore( 'glGetUniformIndices', program, uniformCount, uint8(const), uniformIndices );

return
