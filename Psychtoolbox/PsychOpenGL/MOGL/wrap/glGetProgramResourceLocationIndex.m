function r = glGetProgramResourceLocationIndex( program, programInterface, name )

% glGetProgramResourceLocationIndex  Interface to OpenGL function glGetProgramResourceLocationIndex
%
% usage:  r = glGetProgramResourceLocationIndex( program, programInterface, name )
%
% C function:  GLint glGetProgramResourceLocationIndex(GLuint program, GLenum programInterface, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

r = moglcore( 'glGetProgramResourceLocationIndex', program, programInterface, uint8(name) );

return
