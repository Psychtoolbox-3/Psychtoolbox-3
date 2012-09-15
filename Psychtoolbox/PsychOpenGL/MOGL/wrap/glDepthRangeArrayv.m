function glDepthRangeArrayv( first, count, v )

% glDepthRangeArrayv  Interface to OpenGL function glDepthRangeArrayv
%
% usage:  glDepthRangeArrayv( first, count, v )
%
% C function:  void glDepthRangeArrayv(GLuint first, GLsizei count, const GLdouble* v)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=3,
    error('invalid number of arguments');
end

moglcore( 'glDepthRangeArrayv', first, count, double(v) );

return
