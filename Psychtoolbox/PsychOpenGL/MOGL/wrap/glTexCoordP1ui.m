function glTexCoordP1ui( type, coords )

% glTexCoordP1ui  Interface to OpenGL function glTexCoordP1ui
%
% usage:  glTexCoordP1ui( type, coords )
%
% C function:  void glTexCoordP1ui(GLenum type, GLuint coords)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glTexCoordP1ui', type, coords );

return
