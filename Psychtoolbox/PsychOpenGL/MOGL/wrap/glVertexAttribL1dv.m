function glVertexAttribL1dv( index, v )

% glVertexAttribL1dv  Interface to OpenGL function glVertexAttribL1dv
%
% usage:  glVertexAttribL1dv( index, v )
%
% C function:  void glVertexAttribL1dv(GLuint index, const GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL1dv', index, double(v) );

return
