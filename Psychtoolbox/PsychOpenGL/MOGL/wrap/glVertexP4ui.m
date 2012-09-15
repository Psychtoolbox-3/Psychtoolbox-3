function glVertexP4ui( type, value )

% glVertexP4ui  Interface to OpenGL function glVertexP4ui
%
% usage:  glVertexP4ui( type, value )
%
% C function:  void glVertexP4ui(GLenum type, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP4ui', type, value );

return
