function glVertexBindingDivisor( bindingindex, divisor )

% glVertexBindingDivisor  Interface to OpenGL function glVertexBindingDivisor
%
% usage:  glVertexBindingDivisor( bindingindex, divisor )
%
% C function:  void glVertexBindingDivisor(GLuint bindingindex, GLuint divisor)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexBindingDivisor', bindingindex, divisor );

return
