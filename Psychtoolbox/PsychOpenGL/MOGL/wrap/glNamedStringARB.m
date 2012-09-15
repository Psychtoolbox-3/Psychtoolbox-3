function glNamedStringARB( type, namelen, name, stringlen, string )

% glNamedStringARB  Interface to OpenGL function glNamedStringARB
%
% usage:  glNamedStringARB( type, namelen, name, stringlen, string )
%
% C function:  void glNamedStringARB(GLenum type, GLint namelen, const GLchar* name, GLint stringlen, const GLchar* string)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=5,
    error('invalid number of arguments');
end

moglcore( 'glNamedStringARB', type, namelen, uint8(name), stringlen, uint8(string) );

return
