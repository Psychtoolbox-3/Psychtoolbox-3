function params = glGetFramebufferAttachmentParameteriv( target, attachment, pname )

% glGetFramebufferAttachmentParameteriv  Interface to OpenGL function glGetFramebufferAttachmentParameteriv
%
% usage:  params = glGetFramebufferAttachmentParameteriv( target, attachment, pname )
%
% C function:  void glGetFramebufferAttachmentParameteriv(GLenum target, GLenum attachment, GLenum pname, GLint* params)

% 30-Aug-2012 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=3,
    error('invalid number of arguments');
end

params = int32(0);

moglcore( 'glGetFramebufferAttachmentParameteriv', target, attachment, pname, params );

return
