function glColor4f( red, green, blue, alpha )

% glColor4f  Interface to OpenGL function glColor4f
%
% usage:  glColor4f( red, green, blue, alpha )
%
% C function:  void glColor4f(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---protected---

if nargin~=4,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glColor4f', red, green, blue, alpha );
else
    moglcore( 'ftglColor4f', red, green, blue, alpha );
end

return
