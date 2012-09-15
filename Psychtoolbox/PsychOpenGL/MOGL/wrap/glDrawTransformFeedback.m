function glDrawTransformFeedback( mode, id )

% glDrawTransformFeedback  Interface to OpenGL function glDrawTransformFeedback
%
% usage:  glDrawTransformFeedback( mode, id )
%
% C function:  void glDrawTransformFeedback(GLenum mode, GLuint id)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glDrawTransformFeedback', mode, id );

return
