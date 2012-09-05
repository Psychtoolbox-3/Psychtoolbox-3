function r = gluScaleImage( format, wIn, hIn, typeIn, dataIn, wOut, hOut, typeOut, dataOut )

% gluScaleImage  Interface to OpenGL function gluScaleImage
%
% usage:  r = gluScaleImage( format, wIn, hIn, typeIn, dataIn, wOut, hOut, typeOut, dataOut )
%
% C function:  GLint gluScaleImage(GLenum format, GLsizei wIn, GLsizei hIn, GLenum typeIn, const void* dataIn, GLsizei wOut, GLsizei hOut, GLenum typeOut, GLvoid* dataOut)

% 26-Mar-2011 -- created (generated automatically from header files)

% ---allocate---
% ---protected---

if nargin~=9,
    error('invalid number of arguments');
end

r = moglcore( 'gluScaleImage', format, wIn, hIn, typeIn, dataIn, wOut, hOut, typeOut, dataOut );

return
