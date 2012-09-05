function glVertexAttribL3dv( index, v )

% glVertexAttribL3dv  Interface to OpenGL function glVertexAttribL3dv
%
% usage:  glVertexAttribL3dv( index, v )
%
% C function:  void glVertexAttribL3dv(GLuint index, const GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL3dv', index, double(v) );

return
