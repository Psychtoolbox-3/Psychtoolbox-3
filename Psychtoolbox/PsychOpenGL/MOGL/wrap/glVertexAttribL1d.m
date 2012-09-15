function glVertexAttribL1d( index, x )

% glVertexAttribL1d  Interface to OpenGL function glVertexAttribL1d
%
% usage:  glVertexAttribL1d( index, x )
%
% C function:  void glVertexAttribL1d(GLuint index, GLdouble x)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL1d', index, x );

return
