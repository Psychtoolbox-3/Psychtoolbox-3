function glBlendEquationSeparateiARB( buf, modeRGB, modeAlpha )

% glBlendEquationSeparateiARB  Interface to OpenGL function glBlendEquationSeparateiARB
%
% usage:  glBlendEquationSeparateiARB( buf, modeRGB, modeAlpha )
%
% C function:  void glBlendEquationSeparateiARB(GLuint buf, GLenum modeRGB, GLenum modeAlpha)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glBlendEquationSeparateiARB', buf, modeRGB, modeAlpha );

return
