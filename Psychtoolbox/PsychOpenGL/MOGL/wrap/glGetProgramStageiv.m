function values = glGetProgramStageiv( program, shadertype, pname )

% glGetProgramStageiv  Interface to OpenGL function glGetProgramStageiv
%
% usage:  values = glGetProgramStageiv( program, shadertype, pname )
%
% C function:  void glGetProgramStageiv(GLuint program, GLenum shadertype, GLenum pname, GLint* values)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

values = int32(0);

moglcore( 'glGetProgramStageiv', program, shadertype, pname, values );

return
