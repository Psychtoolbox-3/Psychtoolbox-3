function params = glGetActiveUniformBlockiv( program, uniformBlockIndex, pname )

% glGetActiveUniformBlockiv  Interface to OpenGL function glGetActiveUniformBlockiv
%
% usage:  params = glGetActiveUniformBlockiv( program, uniformBlockIndex, pname )
%
% C function:  void glGetActiveUniformBlockiv(GLuint program, GLuint uniformBlockIndex, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

global GL;

if nargin~=3,
    error('invalid number of arguments');
end

% GL.UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES has a variable number of return indices, which
% we need to query to prealloc properly:
if pname == GL.UNIFORM_BLOCK_ACTIVE_UNIFORM_INDICES
    params = int32(zeros(1, glGetActiveUniformBlockiv( program, uniformBlockIndex, GL.UNIFORM_BLOCK_ACTIVE_UNIFORMS)));
else
    params = int32(0);
end

moglcore( 'glGetActiveUniformBlockiv', program, uniformBlockIndex, pname, params );

return
