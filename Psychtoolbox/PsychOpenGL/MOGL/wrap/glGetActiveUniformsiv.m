function params = glGetActiveUniformsiv( program, uniformCount, uniformIndices, pname )

% glGetActiveUniformsiv  Interface to OpenGL function glGetActiveUniformsiv
%
% usage:  params = glGetActiveUniformsiv( program, uniformCount, uniformIndices, pname )
%
% C function:  void glGetActiveUniformsiv(GLuint program, GLsizei uniformCount, const GLuint* uniformIndices, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

params = int32(zeros(1,uniformCount));

moglcore( 'glGetActiveUniformsiv', program, uniformCount, uint32(uniformIndices), pname, params );

return
