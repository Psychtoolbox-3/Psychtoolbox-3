function glVertexAttribL2dv( index, v )

% glVertexAttribL2dv  Interface to OpenGL function glVertexAttribL2dv
%
% usage:  glVertexAttribL2dv( index, v )
%
% C function:  void glVertexAttribL2dv(GLuint index, const GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL2dv', index, double(v) );

return
