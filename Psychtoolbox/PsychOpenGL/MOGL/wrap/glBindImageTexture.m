function glBindImageTexture( unit, texture, level, layered, layer, access, format )

% glBindImageTexture  Interface to OpenGL function glBindImageTexture
%
% usage:  glBindImageTexture( unit, texture, level, layered, layer, access, format )
%
% C function:  void glBindImageTexture(GLuint unit, GLuint texture, GLint level, GLboolean layered, GLint layer, GLenum access, GLenum format)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=7,
    error('invalid number of arguments');
end

moglcore( 'glBindImageTexture', unit, texture, level, layered, layer, access, format );

return
