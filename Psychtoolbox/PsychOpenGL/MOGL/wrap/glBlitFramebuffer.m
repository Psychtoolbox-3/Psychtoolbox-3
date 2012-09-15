function glBlitFramebuffer( srcX0, srcY0, srcX1, srcY1, dstX0, dstY0, dstX1, dstY1, mask, filter )

% glBlitFramebuffer  Interface to OpenGL function glBlitFramebuffer
%
% usage:  glBlitFramebuffer( srcX0, srcY0, srcX1, srcY1, dstX0, dstY0, dstX1, dstY1, mask, filter )
%
% C function:  void glBlitFramebuffer(GLint srcX0, GLint srcY0, GLint srcX1, GLint srcY1, GLint dstX0, GLint dstY0, GLint dstX1, GLint dstY1, GLbitfield mask, GLenum filter)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=10,
    error('invalid number of arguments');
end

moglcore( 'glBlitFramebuffer', srcX0, srcY0, srcX1, srcY1, dstX0, dstY0, dstX1, dstY1, mask, filter );

return
