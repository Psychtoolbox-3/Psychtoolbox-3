function r = glGetSubroutineIndex( program, shadertype, name )

% glGetSubroutineIndex  Interface to OpenGL function glGetSubroutineIndex
%
% usage:  r = glGetSubroutineIndex( program, shadertype, name )
%
% C function:  GLuint glGetSubroutineIndex(GLuint program, GLenum shadertype, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

r = moglcore( 'glGetSubroutineIndex', program, shadertype, uint8(name) );

return
