function glVertexAttribL4d( index, x, y, z, w )

% glVertexAttribL4d  Interface to OpenGL function glVertexAttribL4d
%
% usage:  glVertexAttribL4d( index, x, y, z, w )
%
% C function:  void glVertexAttribL4d(GLuint index, GLdouble x, GLdouble y, GLdouble z, GLdouble w)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL4d', index, x, y, z, w );

return
