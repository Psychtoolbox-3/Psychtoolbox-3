function glVertexP3ui( type, value )

% glVertexP3ui  Interface to OpenGL function glVertexP3ui
%
% usage:  glVertexP3ui( type, value )
%
% C function:  void glVertexP3ui(GLenum type, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP3ui', type, value );

return
