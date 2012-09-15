function glDeleteVertexArrays( n, arrays )

% glDeleteVertexArrays  Interface to OpenGL function glDeleteVertexArrays
%
% usage:  glDeleteVertexArrays( n, arrays )
%
% C function:  void glDeleteVertexArrays(GLsizei n, const GLuint* arrays)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteVertexArrays', n, uint32(arrays) );

return
