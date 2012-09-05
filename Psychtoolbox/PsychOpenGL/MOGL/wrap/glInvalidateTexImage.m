function glInvalidateTexImage( texture, level )

% glInvalidateTexImage  Interface to OpenGL function glInvalidateTexImage
%
% usage:  glInvalidateTexImage( texture, level )
%
% C function:  void glInvalidateTexImage(GLuint texture, GLint level)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glInvalidateTexImage', texture, level );

return
