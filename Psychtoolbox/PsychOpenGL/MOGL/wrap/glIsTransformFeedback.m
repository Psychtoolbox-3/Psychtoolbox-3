function r = glIsTransformFeedback( id )

% glIsTransformFeedback  Interface to OpenGL function glIsTransformFeedback
%
% usage:  r = glIsTransformFeedback( id )
%
% C function:  GLboolean glIsTransformFeedback(GLuint id)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

r = moglcore( 'glIsTransformFeedback', id );

return
