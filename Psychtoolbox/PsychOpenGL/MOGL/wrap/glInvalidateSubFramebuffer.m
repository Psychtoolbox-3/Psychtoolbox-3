function glInvalidateSubFramebuffer( target, numAttachments, attachments, x, y, width, height )

% glInvalidateSubFramebuffer  Interface to OpenGL function glInvalidateSubFramebuffer
%
% usage:  glInvalidateSubFramebuffer( target, numAttachments, attachments, x, y, width, height )
%
% C function:  void glInvalidateSubFramebuffer(GLenum target, GLsizei numAttachments, const GLenum* attachments, GLint x, GLint y, GLsizei width, GLsizei height)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glInvalidateSubFramebuffer', target, numAttachments, uint32(attachments), x, y, width, height );

return
