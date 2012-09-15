function glProgramParameteri( program, pname, value )

% glProgramParameteri  Interface to OpenGL function glProgramParameteri
%
% usage:  glProgramParameteri( program, pname, value )
%
% C function:  void glProgramParameteri(GLuint program, GLenum pname, GLint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glProgramParameteri', program, pname, value );

return
