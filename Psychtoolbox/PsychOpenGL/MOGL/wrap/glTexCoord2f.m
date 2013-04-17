function glTexCoord2f( s, t )

% glTexCoord2f  Interface to OpenGL function glTexCoord2f
%
% usage:  glTexCoord2f( s, t )
%
% C function:  void glTexCoord2f(GLfloat s, GLfloat t)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---protected---

if nargin~=2,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glTexCoord2f', s, t );
else
    moglcore( 'ftglTexCoord2f', s, t );
end

return
