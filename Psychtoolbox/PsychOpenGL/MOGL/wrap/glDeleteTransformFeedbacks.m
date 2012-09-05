function glDeleteTransformFeedbacks( n, ids )

% glDeleteTransformFeedbacks  Interface to OpenGL function glDeleteTransformFeedbacks
%
% usage:  glDeleteTransformFeedbacks( n, ids )
%
% C function:  void glDeleteTransformFeedbacks(GLsizei n, const GLuint* ids)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDeleteTransformFeedbacks', n, uint32(ids) );

return
