function glSecondaryColorP3ui( type, color )

% glSecondaryColorP3ui  Interface to OpenGL function glSecondaryColorP3ui
%
% usage:  glSecondaryColorP3ui( type, color )
%
% C function:  void glSecondaryColorP3ui(GLenum type, GLuint color)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glSecondaryColorP3ui', type, color );

return
