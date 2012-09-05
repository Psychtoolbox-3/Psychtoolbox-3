function params = glGetProgramInterfaceiv( program, programInterface, pname )

% glGetProgramInterfaceiv  Interface to OpenGL function glGetProgramInterfaceiv
%
% usage:  params = glGetProgramInterfaceiv( program, programInterface, pname )
%
% C function:  void glGetProgramInterfaceiv(GLuint program, GLenum programInterface, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetProgramInterfaceiv', program, programInterface, pname, params );

return
