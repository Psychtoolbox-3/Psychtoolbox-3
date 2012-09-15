function glInvalidateFramebuffer( target, numAttachments, attachments )

% glInvalidateFramebuffer  Interface to OpenGL function glInvalidateFramebuffer
%
% usage:  glInvalidateFramebuffer( target, numAttachments, attachments )
%
% C function:  void glInvalidateFramebuffer(GLenum target, GLsizei numAttachments, const GLenum* attachments)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glInvalidateFramebuffer', target, numAttachments, uint32(attachments) );

return
