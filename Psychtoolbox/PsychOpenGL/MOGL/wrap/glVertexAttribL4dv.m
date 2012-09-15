function glVertexAttribL4dv( index, v )

% glVertexAttribL4dv  Interface to OpenGL function glVertexAttribL4dv
%
% usage:  glVertexAttribL4dv( index, v )
%
% C function:  void glVertexAttribL4dv(GLuint index, const GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexAttribL4dv', index, double(v) );

return
