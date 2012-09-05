function params = glGetProgramiv( program, pname )

% glGetProgramiv  Interface to OpenGL function glGetProgramiv
%
% usage:  params = glGetProgramiv( program, pname )
%
% C function:  void glGetProgramiv(GLuint program, GLenum pname, GLint* params)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetProgramiv', program, pname, params );

return
