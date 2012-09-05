function glVertexAttribL3d( index, x, y, z )

% glVertexAttribL3d  Interface to OpenGL function glVertexAttribL3d
%
% usage:  glVertexAttribL3d( index, x, y, z )
%
% C function:  void glVertexAttribL3d(GLuint index, GLdouble x, GLdouble y, GLdouble z)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=4,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL3d', index, x, y, z );

return
