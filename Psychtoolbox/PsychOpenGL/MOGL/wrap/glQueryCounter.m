function glQueryCounter( id, target )

% glQueryCounter  Interface to OpenGL function glQueryCounter
%
% usage:  glQueryCounter( id, target )
%
% C function:  void glQueryCounter(GLuint id, GLenum target)

% 30-Aug-2012 -- created (generated automatically from header files)

if nargin~=2,
    error('invalid number of arguments');
end

moglcore( 'glQueryCounter', id, target );

return
