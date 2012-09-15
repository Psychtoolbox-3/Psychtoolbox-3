function glVertexAttribDivisorARB( index, divisor )

% glVertexAttribDivisorARB  Interface to OpenGL function glVertexAttribDivisorARB
%
% usage:  glVertexAttribDivisorARB( index, divisor )
%
% C function:  void glVertexAttribDivisorARB(GLuint index, GLuint divisor)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribDivisorARB', index, divisor );

return
