function glEnd

% glEnd  Interface to OpenGL function glEnd
%
% usage:  glEnd
%
% C function:  void glEnd(void)

% 25-Mar-2011 -- created (generated automatically from header files)

% ---protected---

if nargin~=0,
    error('invalid number of arguments');
end

if ~IsGLES
    moglcore( 'glEnd' );
else
    moglcore( 'ftglEnd' );
end

return
