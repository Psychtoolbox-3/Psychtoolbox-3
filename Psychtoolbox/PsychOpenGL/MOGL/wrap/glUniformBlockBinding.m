function glUniformBlockBinding( program, uniformBlockIndex, uniformBlockBinding )

% glUniformBlockBinding  Interface to OpenGL function glUniformBlockBinding
%
% usage:  glUniformBlockBinding( program, uniformBlockIndex, uniformBlockBinding )
%
% C function:  void glUniformBlockBinding(GLuint program, GLuint uniformBlockIndex, GLuint uniformBlockBinding)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glUniformBlockBinding', program, uniformBlockIndex, uniformBlockBinding );

return
