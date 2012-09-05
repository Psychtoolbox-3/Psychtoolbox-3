function glDrawTransformFeedbackStream( mode, id, stream )

% glDrawTransformFeedbackStream  Interface to OpenGL function glDrawTransformFeedbackStream
%
% usage:  glDrawTransformFeedbackStream( mode, id, stream )
%
% C function:  void glDrawTransformFeedbackStream(GLenum mode, GLuint id, GLuint stream)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDrawTransformFeedbackStream', mode, id, stream );

return
