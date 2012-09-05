function glVertexAttribL2d( index, x, y )

% glVertexAttribL2d  Interface to OpenGL function glVertexAttribL2d
%
% usage:  glVertexAttribL2d( index, x, y )
%
% C function:  void glVertexAttribL2d(GLuint index, GLdouble x, GLdouble y)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL2d', index, x, y );

return
