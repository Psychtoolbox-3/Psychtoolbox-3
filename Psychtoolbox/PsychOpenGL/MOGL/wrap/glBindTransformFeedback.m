function glBindTransformFeedback( target, id )

% glBindTransformFeedback  Interface to OpenGL function glBindTransformFeedback
%
% usage:  glBindTransformFeedback( target, id )
%
% C function:  void glBindTransformFeedback(GLenum target, GLuint id)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glBindTransformFeedback', target, id );

return
