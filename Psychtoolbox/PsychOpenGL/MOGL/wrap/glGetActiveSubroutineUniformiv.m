function values = glGetActiveSubroutineUniformiv( program, shadertype, index, pname )

% glGetActiveSubroutineUniformiv  Interface to OpenGL function glGetActiveSubroutineUniformiv
%
% usage:  values = glGetActiveSubroutineUniformiv( program, shadertype, index, pname )
%
% C function:  void glGetActiveSubroutineUniformiv(GLuint program, GLenum shadertype, GLuint index, GLenum pname, GLint* values)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

global GL;

if nargin~=4,
    error('invalid number of arguments');
end

% GL.COMPATIBLE_SUBROUTINES has a variable number of return indices, which
% we need to query to prealloc properly:
if pname == GL.COMPATIBLE_SUBROUTINES
    values = int32(zeros(1, glGetActiveSubroutineUniformiv( program, shadertype, index, GL.NUM_COMPATIBLE_SUBROUTINES)));
else
    values = int32(0);
end

moglcore( 'glGetActiveSubroutineUniformiv', program, shadertype, index, pname, values );

return
