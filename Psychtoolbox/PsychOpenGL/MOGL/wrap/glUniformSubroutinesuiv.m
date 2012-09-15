function glUniformSubroutinesuiv( shadertype, count, indices )

% glUniformSubroutinesuiv  Interface to OpenGL function glUniformSubroutinesuiv
%
% usage:  glUniformSubroutinesuiv( shadertype, count, indices )
%
% C function:  void glUniformSubroutinesuiv(GLenum shadertype, GLsizei count, const GLuint* indices)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUniformSubroutinesuiv', shadertype, count, uint32(indices) );

return
