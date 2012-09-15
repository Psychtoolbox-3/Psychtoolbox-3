function glVertexP2ui( type, value )

% glVertexP2ui  Interface to OpenGL function glVertexP2ui
%
% usage:  glVertexP2ui( type, value )
%
% C function:  void glVertexP2ui(GLenum type, GLuint value)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glVertexP2ui', type, value );

return
