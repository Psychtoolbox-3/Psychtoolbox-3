function glMinSampleShadingARB( value )

% glMinSampleShadingARB  Interface to OpenGL function glMinSampleShadingARB
%
% usage:  glMinSampleShadingARB( value )
%
% C function:  void glMinSampleShadingARB(GLfloat value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=1,
    error('invalid number of arguments');
end

moglcore( 'glMinSampleShadingARB', value );

return
