function glColorP3ui( type, color )

% glColorP3ui  Interface to OpenGL function glColorP3ui
%
% usage:  glColorP3ui( type, color )
%
% C function:  void glColorP3ui(GLenum type, GLuint color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glColorP3ui', type, color );

return
