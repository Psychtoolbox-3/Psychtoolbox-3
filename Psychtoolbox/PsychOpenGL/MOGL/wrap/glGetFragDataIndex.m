function r = glGetFragDataIndex( program, name )

% glGetFragDataIndex  Interface to OpenGL function glGetFragDataIndex
%
% usage:  r = glGetFragDataIndex( program, name )
%
% C function:  GLint glGetFragDataIndex(GLuint program, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

r = moglcore( 'glGetFragDataIndex', program, uint8(name) );

return
