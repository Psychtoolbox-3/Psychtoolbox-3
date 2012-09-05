function glColorP4ui( type, color )

% glColorP4ui  Interface to OpenGL function glColorP4ui
%
% usage:  glColorP4ui( type, color )
%
% C function:  void glColorP4ui(GLenum type, GLuint color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glColorP4ui', type, color );

return
