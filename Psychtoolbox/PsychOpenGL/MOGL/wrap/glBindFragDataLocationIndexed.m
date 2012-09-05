function glBindFragDataLocationIndexed( program, colorNumber, index, name )

% glBindFragDataLocationIndexed  Interface to OpenGL function glBindFragDataLocationIndexed
%
% usage:  glBindFragDataLocationIndexed( program, colorNumber, index, name )
%
% C function:  void glBindFragDataLocationIndexed(GLuint program, GLuint colorNumber, GLuint index, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glBindFragDataLocationIndexed', program, colorNumber, index, uint8(name) );

return
