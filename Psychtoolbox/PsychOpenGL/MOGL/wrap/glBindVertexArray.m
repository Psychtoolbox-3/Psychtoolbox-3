function glBindVertexArray( array )

% glBindVertexArray  Interface to OpenGL function glBindVertexArray
%
% usage:  glBindVertexArray( array )
%
% C function:  void glBindVertexArray(GLuint array)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glBindVertexArray', array );

return
