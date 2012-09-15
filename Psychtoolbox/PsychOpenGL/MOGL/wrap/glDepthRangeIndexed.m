function glDepthRangeIndexed( index, n, f )

% glDepthRangeIndexed  Interface to OpenGL function glDepthRangeIndexed
%
% usage:  glDepthRangeIndexed( index, n, f )
%
% C function:  void glDepthRangeIndexed(GLuint index, GLdouble n, GLdouble f)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDepthRangeIndexed', index, n, f );

return
