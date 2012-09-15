function glTexImage2D( target, level, internalformat, width, height, border, format, type, pixels )

% glTexImage2D  Interface to OpenGL function glTexImage2D
%
% usage:  glTexImage2D( target, level, internalformat, width, height, border, format, type, pixels )
%
% C function:  void glTexImage2D(GLenum target, GLint level, GLenum internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid* pixels)

% 05-Mar-2006 -- created (generated automatically from header files)

% ---protected---

if nargin~=9,
    error('invalid number of arguments');
end

% Special case NULL-Ptr passed?
% if pixels == 0
%     % We have to allocate a matching dummy matrix for now, until
%     % we implement a proper NULL-Ptr interface in moglcore for this
%     % special case - Ugly ugly ugly...:
%     pixels = zeros(width, height, 2);
% end;

moglcore( 'glTexImage2D', target, level, internalformat, width, height, border, format, type, pixels );

return
% ---skip---
