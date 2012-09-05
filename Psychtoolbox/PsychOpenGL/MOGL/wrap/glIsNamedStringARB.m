function r = glIsNamedStringARB( namelen, name )

% glIsNamedStringARB  Interface to OpenGL function glIsNamedStringARB
%
% usage:  r = glIsNamedStringARB( namelen, name )
%
% C function:  GLboolean glIsNamedStringARB(GLint namelen, const GLchar* name)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

r = moglcore( 'glIsNamedStringARB', namelen, uint8(name) );

return
