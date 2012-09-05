function r = glIsRenderbuffer( renderbuffer )

% glIsRenderbuffer  Interface to OpenGL function glIsRenderbuffer
%
% usage:  r = glIsRenderbuffer( renderbuffer )
%
% C function:  GLboolean glIsRenderbuffer(GLuint renderbuffer)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

r = moglcore( 'glIsRenderbuffer', renderbuffer );

return
