function r = glGetProgramResourceLocation( program, programInterface, name )

% glGetProgramResourceLocation  Interface to OpenGL function glGetProgramResourceLocation
%
% usage:  r = glGetProgramResourceLocation( program, programInterface, name )
%
% C function:  GLint glGetProgramResourceLocation(GLuint program, GLenum programInterface, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

r = moglcore( 'glGetProgramResourceLocation', program, programInterface, uint8(name) );

return
