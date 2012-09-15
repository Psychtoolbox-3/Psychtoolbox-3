function r = glGetProgramResourceIndex( program, programInterface, name )

% glGetProgramResourceIndex  Interface to OpenGL function glGetProgramResourceIndex
%
% usage:  r = glGetProgramResourceIndex( program, programInterface, name )
%
% C function:  GLuint glGetProgramResourceIndex(GLuint program, GLenum programInterface, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

r = moglcore( 'glGetProgramResourceIndex', program, programInterface, uint8(name) );

return
