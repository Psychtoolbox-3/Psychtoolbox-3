function glDrawTransformFeedbackInstanced( mode, id, instancecount )

% glDrawTransformFeedbackInstanced  Interface to OpenGL function glDrawTransformFeedbackInstanced
%
% usage:  glDrawTransformFeedbackInstanced( mode, id, instancecount )
%
% C function:  void glDrawTransformFeedbackInstanced(GLenum mode, GLuint id, GLsizei instancecount)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDrawTransformFeedbackInstanced', mode, id, instancecount );

return
